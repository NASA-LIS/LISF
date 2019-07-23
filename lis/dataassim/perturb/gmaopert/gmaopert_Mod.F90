!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !MODULE: gmaopert_Mod
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
module gmaopert_Mod
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_mpiMod
  use LIS_historyMod
  use LIS_fileIOMod
  use LIS_timeMgrMod
  use LIS_DAobservationsMod
  use landpert_module
  use landpert_routines
  use random_fields
  use LIS_ran2_gasdev

  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: gmaoPert_init
  PUBLIC :: gmaoPert_setup
  PUBLIC :: gmaoperturb
  PUBLIC :: gmaopert_readrestart
  PUBLIC :: gmaopert_writerestart

  type forcpertdec
     integer                                     :: N_x, N_y
     type(pert_param_type), dimension(:),pointer :: forcepert_param
     integer, dimension(:,:,:,:), pointer :: Forcepert_rseed
     real, dimension(:,:,:,:), pointer:: Forcepert
     real, dimension(:,:,:,:), pointer:: Forcepert_ntrmdt
  end type forcpertdec
  
  type progpertdec
     integer                                     :: N_x, N_y
     type(pert_param_type), dimension(:),pointer :: progpert_param
     integer, dimension(:,:,:,:), pointer :: progpert_rseed
     real, dimension(:,:,:,:), pointer:: progpert
     real, dimension(:,:,:,:), pointer:: progpert_ntrmdt
  end type progpertdec

  type obspertdec
     integer                                     :: N_x, N_y
     type(pert_param_type), dimension(:),pointer :: obspert_param
     integer, dimension(:,:,:,:), pointer :: obspert_rseed
     real, dimension(:,:,:,:), pointer:: obspert
     real, dimension(:,:,:,:), pointer:: obspert_ntrmdt
  end type obspertdec

  type(forcpertdec), pointer :: forcPert(:,:)
  type(progpertdec), pointer :: progPert(:,:)
  type(obspertdec), pointer :: obsPert(:,:)  
  integer         , pointer :: nobs(:,:)
  logical                   :: f_xyCorr
  logical                   :: p_xyCorr
  logical                   :: o_xyCorr

  
  contains
!BOP
! !ROUTINE: gmaopert_init
!  \label{gmaopert_init}
!
! !DESCRIPTION:
!  Initializes memory structures for specified perturbation
!  options.
! 
! !INTERFACE:      
    subroutine gmaoPert_init(index)
! !USES: 

! !ARGUMENTS: 
      integer, intent(in) :: index
!EOP
      
      if(index.eq.1) then          
         allocate(forcPert(LIS_rc%nnest, 1))            
         f_xycorr = .false. 
      elseif(index.eq.2) then 
         allocate(progPert(LIS_rc%nnest, LIS_rc%nperts))
         p_xycorr = .false. 
      elseif(index.eq.3) then 
         allocate(obspert(LIS_rc%nnest, LIS_rc%ndas))
         allocate(nobs(LIS_rc%nnest,LIS_rc%ndas))
         o_xycorr = .false. 
      endif

    end subroutine gmaoPert_init

!BOP
! !ROUTINE: gmaopert_setup
!  \label{gmaopert_setup}
!
! !INTERFACE:      
    subroutine gmaoPert_setup(index, k, Base_State, Pert_State)
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
      integer                      :: gid
      integer                      :: status
      integer, parameter           :: N_domain = 1
      character*100                :: sname
      integer                      :: objcount
      integer, allocatable             :: init_Pert_rseed(:,:,:,:)
      integer, allocatable             :: ens_id(:)
      integer                          :: domain_id
      character*100,       allocatable :: pertobjs(:)
      type(ESMF_Field),    allocatable :: pertField(:)
      real,    allocatable             :: std(:,:,:) 
      real,    allocatable             :: std1(:)
      integer                          :: i,t,m,tindex,col,row
      integer                          :: c,r
      integer,  allocatable            :: gindex(:,:)
      integer                          :: n, nest
      real                             :: dx, dy, dtstep
      real,     allocatable            :: global_std(:,:)
      integer                          :: global_flag

      call ESMF_StateGet(Base_State(1),name=sname, rc=status)
      call LIS_verify(status, &
           "ESMF_StateGet failed in gmaopert_setup")

      if(index.eq.1) then            
         do nest=1,LIS_rc%nnest
            if(.not.isXYcorrEnabled(Pert_State(nest))) then 
               allocate(ens_id(LIS_rc%nensem(nest)))
               call LIS_getDomainResolutions(nest,dx,dy)
               dtstep = LIS_rc%ts
               ! allocate substructures only if needed. 
               nullify(forcPert(nest,k)%forcepert_param)
               call assemble_pert_param(Pert_State(nest), &
                    LIS_rc%lnc(nest),LIS_rc%lnr(nest),&
                    LIS_rc%nforcepert,forcPert(nest,k)%forcepert_param)
               
               allocate(forcPert(nest,k)%Forcepert_ntrmdt(LIS_rc%nforcepert,&
                    LIS_rc%lnc(nest),LIS_rc%lnr(nest),LIS_rc%nensem(nest)))
               allocate(forcPert(nest,k)%Forcepert(LIS_rc%nforcepert,&
                    LIS_rc%lnc(nest),LIS_rc%lnr(nest),LIS_rc%nensem(nest)))
               ! initialize just in case (should not be needed)      
               forcPert(nest,k)%Forcepert_ntrmdt = 0.    
               forcPert(nest,k)%Forcepert        = 0.   
               
               do n=1,LIS_rc%nensem(nest)
                  ens_id(n) = n-1
               end do
               
               forcPert(nest,k)%N_x = LIS_rc%lnc(nest)
               forcPert(nest,k)%N_y = LIS_rc%lnr(nest)
               
               allocate(forcPert(nest,k)%Forcepert_rseed(NRANDSEED,&
                    forcPert(nest,k)%N_x, forcPert(nest,k)%N_y,&
                    LIS_rc%nensem(nest)))
               
               allocate(init_Pert_rseed(LIS_rc%nensem(nest),&
                    forcPert(nest,k)%N_x, forcPert(nest,k)%N_y,&
                    N_domain))
               allocate(gindex(forcPert(nest,k)%N_x, forcPert(nest,k)%N_y))
               gindex = -1
               
               do r=1,forcPert(nest,k)%N_y
                  do c=1,forcPert(nest,k)%N_x
                     if(r.le.LIS_rc%lnr(nest)) then 
                        row = r + LIS_nss_halo_ind(nest,LIS_localPet+1)-1
                     else
                        row = r - LIS_rc%lnr(nest)
                     endif
                     if(c.le.LIS_rc%lnc(nest)) then 
                        col = c + LIS_ews_halo_ind(nest,LIS_localPet+1)-1
                     else
                        col = c -  LIS_rc%lnc(nest)
                     endif
                     gindex(c,r) = col + (row-1)*LIS_rc%gnc(nest)
                  enddo
               enddo
               
               ! get different negative integer for each ensemble member and
               ! each domain

               do r=1,forcPert(nest,k)%N_y
                  do c=1,forcPert(nest,k)%N_x
                     domain_id = gindex(c,r)
                     call get_init_Pert_rseed(LIS_rc%nensem(nest), &
                          N_domain, LIS_rc%gnc(nest), &
                          LIS_rc%gnr(nest),  &
                          ens_id, domain_id, &
                          init_Pert_rseed (:,c,r,:))
                  enddo
               enddo

               deallocate(gindex)
               
               ! initialize first row of Forcepert_rseed (for first domain)
               do n=1,LIS_rc%nensem(nest)
                  forcPert(nest,k)%forcepert_rseed(1,:,:,n) = &
                       init_Pert_rseed( n, :,:,1)
               enddo
               ! initial call to get_pert  (for first domain)
               !
               ! this initializes Forcepert_ntrmdt and the rest of Forcepert_rseed 
               !
               ! NOTE: after initial call to get_pert use restart files 
               !       for Forcepert_rseed and Forcepert_ntrmdt to continue
               !       the perturbation time series whenever the land model
               !       integration is interrupted
               
               call ESMF_StateGet(Pert_State(nest), itemCount=objcount,rc=status)
               call LIS_verify(status, &
                    "ESMF_StateGet failed in gmaopert_setup")
               
               allocate(pertobjs(objcount))
               allocate(pertField(objcount))
               
               call ESMF_StateGet(Pert_State(nest), itemNameList=pertobjs,&
                    rc=status)
               call LIS_verify(status,&
                    "ESMF_StateGet failed in gmaopert_setup")
               
               allocate(std(objcount,LIS_rc%lnc(nest), LIS_rc%lnr(nest)))
               allocate(std1(LIS_rc%ngrid(nest)))
               
               do i=1,objcount
                  call ESMF_StateGet(Pert_State(nest), pertobjs(i),pertField(i),&
                       rc=status)
                  call LIS_verify(status)
                  
                  if(LIS_rc%ngrid(nest).gt.0) then 
                     call ESMF_AttributeGet(pertField(i),"Standard Deviation",&
                          std1,rc=status)
                     !                       itemCount=LIS_rc%ngrid(nest),rc=status)
                     call LIS_verify(status)
                  endif
                  
                  std(i,:,:) = 0.0 
                  do t=1,LIS_rc%ntiles(nest)/LIS_rc%nensem(nest)            
                     do m=1, LIS_rc%nensem(nest)
                        tindex = (t-1)*LIS_rc%nensem(nest)+m
                        
                        col = LIS_domain(nest)%tile(tindex)%col
                        row = LIS_domain(nest)%tile(tindex)%row
                        gid = LIS_domain(nest)%gindex(col,row)
                        
                        std(i,col,row) = std(i,col,row) + std1(gid)/&
                             LIS_rc%nensem(nest) 
                     enddo
                  enddo
               enddo
               
               call get_pert(                                           &
                    LIS_rc%nforcepert, LIS_rc%nensem(nest),&
                    LIS_rc%lnc(nest),LIS_rc%lnr(nest),  &
                    forcPert(nest,k)%N_x, forcPert(nest,k)%N_y, &
                    dx, dy, dtstep,                                     &
                    std,                              & 
                    forcPert(nest,k)%forcepert_param, &
                    forcPert(nest,k)%Forcepert_rseed, &
                    forcPert(nest,k)%Forcepert_ntrmdt, &
                    forcPert(nest,k)%Forcepert,   &                
                    initialize_rseed=.true.,    &
                    initialize_ntrmdt=.true.            )
               
               deallocate(init_Pert_rseed)
               deallocate(ens_id)
               deallocate(std)
               deallocate(std1)
               deallocate(pertobjs)
               deallocate(pertField)
            else !xy corrs enabled
               f_xyCorr = .true. 

               call LIS_getDomainResolutions(nest,dx,dy)
               dtstep = LIS_rc%ts
               
               nullify(forcPert(nest,k)%forcepert_param)
               call assemble_pert_param(Pert_State(nest), &
                    LIS_rc%gnc(nest),LIS_rc%gnr(nest),&
                    LIS_rc%nforcepert,forcPert(nest,k)%forcepert_param,&
                    global_flag)
               
               if(LIS_masterproc) then
                  allocate(ens_id(LIS_rc%nensem(nest)))              
                  allocate(forcPert(nest,k)%forcepert_ntrmdt(LIS_rc%nforcepert,&
                       LIS_rc%gnc(nest),LIS_rc%gnr(nest),LIS_rc%nensem(nest)))
                  allocate(forcPert(nest,k)%forcepert(LIS_rc%nforcepert,&
                       LIS_rc%gnc(nest),&
                       LIS_rc%gnr(nest),LIS_rc%nensem(nest)))
                  ! initialize just in case (should not be needed)      
                  forcPert(nest,k)%forcepert_ntrmdt = 0.    
                  forcPert(nest,k)%forcepert        = 0.   
                  
                  do n=1,LIS_rc%nensem(nest)
                     ens_id(n) = n-1
                  end do
               
                  call get_fft_grid( LIS_rc%gnc(nest), LIS_rc%gnr(nest), &
                       dx, dy, &
                       forcPert(nest,k)%forcepert_param(1)%xcorr, &
                       forcPert(nest,k)%forcepert_param(1)%ycorr, &
                       forcPert(nest,k)%N_x, forcPert(nest,k)%N_y)
                  
                  allocate(forcPert(nest,k)%forcepert_rseed(NRANDSEED,&
                       forcPert(nest,k)%N_x, forcPert(nest,k)%N_y,&
                       LIS_rc%nensem(nest)))
               
                  allocate(init_Pert_rseed(LIS_rc%nensem(nest),&
                       forcPert(nest,k)%N_x, forcPert(nest,k)%N_y,&
                       N_domain))
                  allocate(gindex(forcPert(nest,k)%N_x, forcPert(nest,k)%N_y))
                  gindex = -1
               
                  do r=1,forcPert(nest,k)%N_y
                     do c=1,forcPert(nest,k)%N_x
                        row = r 
                        col = c 
                        gindex(c,r) = col + (row-1)*LIS_rc%gnc(nest)                        
                     enddo
                  enddo
                  
               
                  ! get different negative integer for each ensemble member and
                  ! each domain
                  
                  do r=1,forcPert(nest,k)%N_y
                     do c=1,forcPert(nest,k)%N_x
                        domain_id = gindex(c,r)
                        call get_init_Pert_rseed(LIS_rc%nensem(nest), &
                             N_domain, LIS_rc%gnc(nest), &
                             LIS_rc%gnr(nest), &
                             ens_id, domain_id, &
                             init_Pert_rseed (:,c,r,:))
                     enddo
                  enddo
                  deallocate(gindex)
               
               ! initialize first row of Forcepert_rseed (for first domain)
                  do n=1,LIS_rc%nensem(nest)
                     forcPert(nest,k)%forcepert_rseed(1,:,:,n) = &
                          init_Pert_rseed( n, :,:,1)
                  enddo
               else
                  allocate(forcPert(nest,k)%forcepert(LIS_rc%nforcepert,&
                       1, 1, LIS_rc%nensem(nest)))
               endif
               ! initial call to get_pert  (for first domain)
               !
               call ESMF_StateGet(Pert_State(nest), itemCount=objcount,rc=status)
               call LIS_verify(status, 'ESMF_StateGet failed in gmaopert_setup')
               
               allocate(pertobjs(objcount))
               allocate(pertField(objcount))
               
               call ESMF_StateGet(Pert_State(nest), itemNameList=pertobjs, rc=status)
               call LIS_verify(status, 'ESMF_StateGet failed in gmaopert_setup')
               
               if(LIS_masterproc) then 
                  allocate(std(objcount,LIS_rc%gnc(nest), LIS_rc%gnr(nest)))
                  allocate(std1(LIS_rc%ngrid(nest)))
               else
                  allocate(std(1,1,1))
                  allocate(std1(LIS_rc%ngrid(nest)))
               endif

               do i=1,objcount
                  call ESMF_StateGet(Pert_State(nest), pertobjs(i),pertField(i),&
                       rc=status)
                  call LIS_verify(status, 'ESMF_StateGet failed in gmaopert_setup')
                  
                  if(LIS_rc%ngrid(nest).gt.0) then 
                     call ESMF_AttributeGet(pertField(i),"Standard Deviation",&
                          std1,rc=status)
                     call LIS_verify(status,'ESMF_AttributeGet failed in gmaopert_setup')
                  endif
                  
                  call LIS_gather_1dgrid_to_2dgrid(nest, global_std, std1)
                  
                  if(LIS_masterproc) then 
                     std(i,:,:) = global_std(:,:)
                  endif

                  deallocate(global_std)

               enddo
               
               if(LIS_masterproc) then 
                  call get_pert(                                           &
                       LIS_rc%nforcepert, LIS_rc%nensem(nest),&
                       LIS_rc%gnc(nest),LIS_rc%gnr(nest),  &
                       forcPert(nest,k)%N_x, forcPert(nest,k)%N_y, &
                       dx, dy, dtstep,                                     &
                       std,                              & 
                       forcPert(nest,k)%forcepert_param, &
                       forcPert(nest,k)%forcepert_rseed, &
                       forcPert(nest,k)%forcepert_ntrmdt, &
                       forcPert(nest,k)%forcepert,   &                
                       initialize_rseed=.true.,    &
                       initialize_ntrmdt=.true.            )
                  
                  deallocate(init_Pert_rseed)
                  deallocate(ens_id)
               endif
               deallocate(std)
               deallocate(std1)
            endif

         enddo
      elseif(index.eq.2) then 
         
         do nest=1,LIS_rc%nnest
            if(.not.isXYcorrEnabled(Pert_State(nest))) then 
               allocate(ens_id(LIS_rc%nensem(nest)))
               
               call LIS_getDomainResolutions(nest,dx,dy)
               dtstep = LIS_rc%ts
               
               nullify(progPert(nest,k)%progpert_param)
               call assemble_pert_param(Pert_State(nest), &
                    LIS_rc%lnc(nest),LIS_rc%lnr(nest),&
                    LIS_rc%nstvars(k),progPert(nest,k)%progpert_param)
               
               allocate(progPert(nest,k)%progpert_ntrmdt(LIS_rc%nstvars(k),&
                    LIS_rc%lnc(nest),LIS_rc%lnr(nest),LIS_rc%nensem(nest)))
               allocate(progPert(nest,k)%Progpert(LIS_rc%nstvars(k),&
                    LIS_rc%lnc(nest),&
                    LIS_rc%lnr(nest),LIS_rc%nensem(nest)))
               ! initialize just in case (should not be needed)      
               progPert(nest,k)%progpert_ntrmdt = 0.    
               progPert(nest,k)%progpert        = 0.   
               
               do n=1,LIS_rc%nensem(nest)
                  ens_id(n) = n-1
               end do
               
               progPert(nest,k)%N_x = LIS_rc%lnc(nest)
               progPert(nest,k)%N_y = LIS_rc%lnr(nest)
               
               allocate(progPert(nest,k)%progpert_rseed(NRANDSEED,&
                    progPert(nest,k)%N_x, progPert(nest,k)%N_y,LIS_rc%nensem(nest)))
               
               allocate(init_Pert_rseed(LIS_rc%nensem(nest),&
                    progPert(nest,k)%N_x, progPert(nest,k)%N_y,&
                    N_domain))
               allocate(gindex(progPert(nest,k)%N_x, progPert(nest,k)%N_y))
               gindex = -1
               
               do r=1,progPert(nest,k)%N_y
                  do c=1,progPert(nest,k)%N_x
                     if(r.le.LIS_rc%lnr(nest)) then 
                        row = r + LIS_nss_halo_ind(nest,LIS_localPet+1)-1
                     else
                        row = r - LIS_rc%lnr(nest)
                     endif
                     if(c.le.LIS_rc%lnc(nest)) then 
                        col = c + LIS_ews_halo_ind(nest,LIS_localPet+1)-1
                     else
                        col = c -  LIS_rc%lnc(nest)
                     endif
                     gindex(c,r) = col + (row-1)*LIS_rc%gnc(nest)
                  enddo
               enddo
               
               
            ! get different negative integer for each ensemble member and
            ! each domain
                  
               do r=1,progPert(nest,k)%N_y
                  do c=1,progPert(nest,k)%N_x
                     domain_id = gindex(c,r)
                     call get_init_Pert_rseed(LIS_rc%nensem(nest), &
                          N_domain, LIS_rc%gnc(nest), &
                          LIS_rc%gnr(nest), &
                          ens_id, domain_id, &
                          init_Pert_rseed (:,c,r,:))
                  enddo
               enddo
               deallocate(gindex)
               
               ! initialize first row of Forcepert_rseed (for first domain)
               do n=1,LIS_rc%nensem(nest)
                  progPert(nest,k)%progpert_rseed(1,:,:,n) = &
                       init_Pert_rseed( n, :,:,1)
               enddo
               
               ! initial call to get_pert  (for first domain)
               !
               ! this initializes progpert_ntrmdt and the rest of progpert_rseed 
               !
               ! NOTE: after initial call to get_pert use restart files 
               !       for progpert_rseed and progpert_ntrmdt to continue
               !       the perturbation time series whenever the land model
               !       integration is interrupted
               
               call ESMF_StateGet(Pert_State(nest), itemCount=objcount,rc=status)
               call LIS_verify(status)
               
               allocate(pertobjs(objcount))
               allocate(pertField(objcount))
               
               call ESMF_StateGet(Pert_State(nest), itemNameList=pertobjs, rc=status)
               call LIS_verify(status)
               
               allocate(std(objcount,LIS_rc%lnc(nest), LIS_rc%lnr(nest)))
               allocate(std1(LIS_rc%ngrid(nest)))
               
               do i=1,objcount
                  call ESMF_StateGet(Pert_State(nest), pertobjs(i),pertField(i),&
                       rc=status)
                  call LIS_verify(status)
                  
                  if(LIS_rc%ngrid(nest).gt.0) then 
                     call ESMF_AttributeGet(pertField(i),"Standard Deviation",&
                          std1,rc=status)
                     !                       itemCount=LIS_rc%ngrid(nest),rc=status)
                     call LIS_verify(status)
                  endif
                  
                  std(i,:,:) = 0.0 
                  do t=1,LIS_rc%npatch(nest,LIS_rc%lsm_index)/LIS_rc%nensem(nest)            
                     do m=1, LIS_rc%nensem(nest)
                        tindex = (t-1)*LIS_rc%nensem(nest)+m
                        
                        col = LIS_surface(nest,1)%tile(tindex)%col
                        row = LIS_surface(nest,1)%tile(tindex)%row
                        gid = LIS_domain(nest)%gindex(col,row)
                        
                        std(i,col,row) = std(i,col,row) + std1(gid)/&
                             LIS_rc%nensem(nest) 
                     enddo
                  enddo
               enddo
               
               call get_pert(                                           &
                    LIS_rc%nstvars(k), LIS_rc%nensem(nest),&
                    LIS_rc%lnc(nest),LIS_rc%lnr(nest),  &
                    progPert(nest,k)%N_x, progPert(nest,k)%N_y, &
                    dx, dy, dtstep,                                     &
                    std,                              & 
                    progPert(nest,k)%progpert_param, &
                    progPert(nest,k)%Progpert_rseed, &
                    progPert(nest,k)%Progpert_ntrmdt, &
                    progPert(nest,k)%Progpert,   &                
                    initialize_rseed=.true.,    &
                    initialize_ntrmdt=.true.            )
               deallocate(init_Pert_rseed)
               deallocate(ens_id)
               deallocate(std)
               deallocate(std1)
            else ! xy corrs enabled
               p_xyCorr = .true. 
               call LIS_getDomainResolutions(nest,dx,dy)
               dtstep = LIS_rc%ts
               
               nullify(progPert(nest,k)%progpert_param)
               call assemble_pert_param(Pert_State(nest), &
                    LIS_rc%gnc(nest),LIS_rc%gnr(nest),&
                    LIS_rc%nstvars(k),progPert(nest,k)%progpert_param,&
                    global_flag)
               
               if(LIS_masterproc) then
                  allocate(ens_id(LIS_rc%nensem(nest)))              
                  allocate(progPert(nest,k)%progpert_ntrmdt(LIS_rc%nstvars(k),&
                       LIS_rc%gnc(nest),LIS_rc%gnr(nest),LIS_rc%nensem(nest)))
                  allocate(progPert(nest,k)%Progpert(LIS_rc%nstvars(k),&
                       LIS_rc%gnc(nest),&
                       LIS_rc%gnr(nest),LIS_rc%nensem(nest)))
               ! initialize just in case (should not be needed)      
                  progPert(nest,k)%progpert_ntrmdt = 0.    
                  progPert(nest,k)%progpert        = 0.   
                  
                  do n=1,LIS_rc%nensem(nest)
                     ens_id(n) = n-1
                  end do
               
                  call get_fft_grid( LIS_rc%gnc(nest), LIS_rc%gnr(nest), &
                       dx, dy, &
                       progPert(nest,k)%progpert_param(1)%xcorr, &
                       progPert(nest,k)%progpert_param(1)%ycorr, &
                       progPert(nest,k)%N_x, progPert(nest,k)%N_y)
                  
                  allocate(progPert(nest,k)%progpert_rseed(NRANDSEED,&
                       progPert(nest,k)%N_x, progPert(nest,k)%N_y,&
                       LIS_rc%nensem(nest)))
               
                  allocate(init_Pert_rseed(LIS_rc%nensem(nest),&
                       progPert(nest,k)%N_x, progPert(nest,k)%N_y,&
                       N_domain))
                  allocate(gindex(progPert(nest,k)%N_x, progPert(nest,k)%N_y))
                  gindex = -1
               
                  do r=1,progPert(nest,k)%N_y
                     do c=1,progPert(nest,k)%N_x
                        row = r 
                        col = c 
                        gindex(c,r) = col + (row-1)*LIS_rc%gnc(nest)                        
                     enddo
                  enddo
                  
               
                  ! get different negative integer for each ensemble member and
                  ! each domain
                  
                  do r=1,progPert(nest,k)%N_y
                     do c=1,progPert(nest,k)%N_x
                        domain_id = gindex(c,r)
                        call get_init_Pert_rseed(LIS_rc%nensem(nest), &
                             N_domain, LIS_rc%gnc(nest), &
                             LIS_rc%gnr(nest), &
                             ens_id, domain_id, &
                             init_Pert_rseed (:,c,r,:))
                     enddo
                  enddo
                  deallocate(gindex)
               
               ! initialize first row of Forcepert_rseed (for first domain)
                  do n=1,LIS_rc%nensem(nest)
                     progPert(nest,k)%progpert_rseed(1,:,:,n) = &
                          init_Pert_rseed( n, :,:,1)
                  enddo
               else
                  allocate(progPert(nest,k)%Progpert(LIS_rc%nstvars(k),&
                       1, 1, LIS_rc%nensem(nest)))
               endif
               ! initial call to get_pert  (for first domain)
               !
               ! this initializes progpert_ntrmdt and the rest of progpert_rseed 
               !
               ! NOTE: after initial call to get_pert use restart files 
               !       for progpert_rseed and progpert_ntrmdt to continue
               !       the perturbation time series whenever the land model
               !       integration is interrupted
               
               call ESMF_StateGet(Pert_State(nest), itemCount=objcount,rc=status)
               call LIS_verify(status)
               
               allocate(pertobjs(objcount))
               allocate(pertField(objcount))
               
               call ESMF_StateGet(Pert_State(nest), itemNameList=pertobjs, rc=status)
               call LIS_verify(status)
               
               if(LIS_masterproc) then 
                  allocate(std(objcount,LIS_rc%gnc(nest), LIS_rc%gnr(nest)))
                  allocate(std1(LIS_rc%ngrid(nest)))
               else
                  allocate(std(1,1,1))
                  allocate(std1(LIS_rc%ngrid(nest)))
               endif

               do i=1,objcount
                  call ESMF_StateGet(Pert_State(nest), pertobjs(i),pertField(i),&
                       rc=status)
                  call LIS_verify(status)
                  
                  if(LIS_rc%ngrid(nest).gt.0) then 
                     call ESMF_AttributeGet(pertField(i),"Standard Deviation",&
                          std1,rc=status)
                     call LIS_verify(status)
                  endif
                  
                  call LIS_gather_1dgrid_to_2dgrid(nest, global_std, std1)
                  
                  if(LIS_masterproc) then 
                     std(i,:,:) = global_std(:,:)
                  endif

                  deallocate(global_std)

               enddo
               
               if(LIS_masterproc) then 
                  call get_pert(                                           &
                       LIS_rc%nstvars(k), LIS_rc%nensem(nest),&
                       LIS_rc%gnc(nest),LIS_rc%gnr(nest),  &
                       progPert(nest,k)%N_x, progPert(nest,k)%N_y, &
                       dx, dy, dtstep,                                     &
                       std,                              & 
                       progPert(nest,k)%progpert_param, &
                       progPert(nest,k)%Progpert_rseed, &
                       progPert(nest,k)%Progpert_ntrmdt, &
                       progPert(nest,k)%Progpert,   &                
                       initialize_rseed=.true.,    &
                       initialize_ntrmdt=.true.            )
                  
                  deallocate(init_Pert_rseed)
                  deallocate(ens_id)
               endif
               deallocate(std)
               deallocate(std1)
            endif
         enddo
      elseif(index.eq.3) then 
         
         do nest=1,LIS_rc%nnest
            if(.not.isXYcorrEnabled(Pert_State(nest))) then 

               allocate(ens_id(LIS_rc%nensem(nest)))
               call ESMF_StateGet(Pert_State(nest), itemCount=objcount,rc=status)
               call LIS_verify(status, 'ESMF_StateGet failed in gmaopert_mod')
               
               nobs(nest,k) = objcount
               
               call LIS_getObsDomainResolutions(k,dx,dy)
               dtstep = LIS_rc%ts
               
               nullify(obspert(nest,k)%obspert_param)
               call assemble_pert_param_obs(k, Pert_State(nest), &
                    LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k),&
                    objcount,obspert(nest,k)%obspert_param)
               
               allocate(obspert(nest,k)%obspert_ntrmdt(objcount,&
                    LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k),LIS_rc%nensem(nest)))
               allocate(obspert(nest,k)%Obspert(objcount,LIS_rc%obs_lnc(k),&
                    LIS_rc%obs_lnr(k),LIS_rc%nensem(nest)))
               ! initialize just in case (should not be needed)      
               obspert(nest,k)%obspert_ntrmdt = 0.    
               obspert(nest,k)%obspert        = 0.   
               
               do n=1,LIS_rc%nensem(nest)
                  ens_id(n) = n-1
               end do
               
               obsPert(nest,k)%N_x = LIS_rc%obs_lnc(k)
               obsPert(nest,k)%N_y = LIS_rc%obs_lnr(k)               
               
               allocate(obspert(nest,k)%obspert_rseed(NRANDSEED,&
                    obsPert(nest,k)%N_x, obsPert(nest,k)%N_y,LIS_rc%nensem(nest)))
               
               allocate(init_Pert_rseed(LIS_rc%nensem(nest),&
                    obsPert(nest,k)%N_x, obsPert(nest,k)%N_y,&
                    N_domain))
               allocate(gindex(obsPert(nest,k)%N_x, obsPert(nest,k)%N_y))
               gindex = -1
               do r=1,obsPert(nest,k)%N_y
                  do c=1,obsPert(nest,k)%N_x
                     if(r.le.LIS_rc%obs_lnr(k)) then 
                        row = r + LIS_nss_obs_halo_ind(k,LIS_localPet+1)-1
                     else
                        row = r - LIS_rc%obs_lnr(k)
                     endif
                     if(c.le.LIS_rc%obs_lnc(k)) then 
                        col = c + LIS_ews_obs_halo_ind(k,LIS_localPet+1)-1
                     else
                        col = c -  LIS_rc%obs_lnc(k)
                     endif

                     gindex(c,r) = col + (row-1)*LIS_rc%obs_gnc(k)
                  enddo
               enddo
               
               ! get different negative integer for each ensemble member and
               ! each domain
               
               do r=1,obsPert(nest,k)%N_y
                  do c=1,obsPert(nest,k)%N_x
                     domain_id = gindex(c,r)
                     call get_init_Pert_rseed(LIS_rc%nensem(nest), &
                          N_domain, LIS_rc%gnc(nest), &
                          LIS_rc%gnr(nest),  &
                          ens_id, domain_id, &
                          init_Pert_rseed (:,c,r,:))
                  enddo
               enddo
               deallocate(gindex)
               
            ! initialize first row of Forcepert_rseed (for first domain)
               do n=1,LIS_rc%nensem(nest)
                  obsPert(nest,k)%obspert_rseed(1,:,:,n) = &
                       init_Pert_rseed( n, :,:,1)
               enddo
               ! initial call to get_pert  (for first domain)
               !
               ! this initializes obspert_ntrmdt and the rest of obspert_rseed 
               !
               ! NOTE: after initial call to get_pert use restart files 
               !       for obspert_rseed and obspert_ntrmdt to continue
               !       the perturbation time series whenever the land model
               !       integration is interrupted
               call ESMF_StateGet(Pert_State(nest), itemCount=objcount,rc=status)
               call LIS_verify(status,'ESMF_StateGet:pert_state1 failed in gmaopert_mod')
               
               allocate(pertobjs(objcount))
               allocate(pertField(objcount))
               call ESMF_StateGet(Pert_State(nest), itemNameList=pertobjs, rc=status)
               call LIS_verify(status,'ESMF_StateGet:pert_state2 failed in gmaopert_mod')
               
               allocate(std(objcount, LIS_rc%obs_lnc(k), LIS_rc%obs_lnr(k)))
               allocate(std1(LIS_rc%obs_ngrid(k)))
               
               do i=1,objcount
                  call ESMF_StateGet(Pert_State(nest), pertobjs(i),pertField(i),&
                       rc=status)
                  call LIS_verify(status,'ESMF_StateGet:pert_state3 failed in gmaopert_mod')
                  
                  if(LIS_rc%obs_ngrid(k).gt.0) then 
                     call ESMF_AttributeGet(pertField(i),"Standard Deviation",&
                          std1,rc=status)
                     call LIS_verify(status,'ESMF_AttributeGet failed in gmaopert_mod')
                  endif
                  
                  std(i,:,:) = 0.0 
                  do t=1,LIS_rc%obs_ngrid(k)
                     col = LIS_obs_domain(nest,k)%col(t)
                     row = LIS_obs_domain(nest,k)%row(t)
                     
                     std(i,col,row) = std1(t)
                  enddo
               enddo
               call get_pert(                                           &
                    objcount, LIS_rc%nensem(nest),&
                    LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k),  &
                    obsPert(nest,k)%N_x, obsPert(nest,k)%N_y, &
                    dx, dy, dtstep,                                     &
                    std,                              & 
                    obspert(nest,k)%obspert_param, &
                    obspert(nest,k)%Obspert_rseed, &
                    obspert(nest,k)%Obspert_ntrmdt, &
                    obspert(nest,k)%Obspert,   &                
                    initialize_rseed=.true.,    &
                    initialize_ntrmdt=.true.            )
               
               deallocate(init_Pert_rseed)
               deallocate(ens_id)
               deallocate(std)
               deallocate(std1)
            else !xy corrs enabled
               o_xyCorr = .true. 
               call LIS_getObsDomainResolutions(k,dx,dy)
               dtstep = LIS_rc%ts

               call ESMF_StateGet(Pert_State(nest), itemCount=objcount,rc=status)
               call LIS_verify(status)
               
               nobs(nest,k) = objcount
               
               nullify(obsPert(nest,k)%obspert_param)
               call assemble_pert_param_obs(k,Pert_State(nest), &
                    LIS_rc%obs_gnc(k),LIS_rc%obs_gnr(k),&
                    objcount,obsPert(nest,k)%obspert_param,&
                    global_flag)
               
               if(LIS_masterproc) then
                  allocate(ens_id(LIS_rc%nensem(nest)))              
                  allocate(obsPert(nest,k)%obspert_ntrmdt(objcount,&
                       LIS_rc%obs_gnc(k),LIS_rc%obs_gnr(k),LIS_rc%nensem(nest)))
                  allocate(obsPert(nest,k)%obspert(objcount,&
                       LIS_rc%obs_gnc(k),&
                       LIS_rc%obs_gnr(k),LIS_rc%nensem(nest)))
               ! initialize just in case (should not be needed)      
                  obsPert(nest,k)%obspert_ntrmdt = 0.    
                  obsPert(nest,k)%obspert        = 0.   
                  
                  do n=1,LIS_rc%nensem(nest)
                     ens_id(n) = n-1
                  end do
               
                  call get_fft_grid( LIS_rc%obs_gnc(k), LIS_rc%obs_gnr(k), &
                       dx, dy, &
                       obsPert(nest,k)%obspert_param(1)%xcorr, &
                       obsPert(nest,k)%obspert_param(1)%ycorr, &
                       obsPert(nest,k)%N_x, obsPert(nest,k)%N_y)
                  
                  allocate(obsPert(nest,k)%obspert_rseed(NRANDSEED,&
                       obsPert(nest,k)%N_x, obsPert(nest,k)%N_y,&
                       LIS_rc%nensem(nest)))
               
                  allocate(init_Pert_rseed(LIS_rc%nensem(nest),&
                       obsPert(nest,k)%N_x, obsPert(nest,k)%N_y,&
                       N_domain))
                  allocate(gindex(obsPert(nest,k)%N_x, obsPert(nest,k)%N_y))
                  gindex = -1
               
                  do r=1,obsPert(nest,k)%N_y
                     do c=1,obsPert(nest,k)%N_x
                        row = r 
                        col = c 
                        gindex(c,r) = col + (row-1)*LIS_rc%obs_gnc(k)                        
                     enddo
                  enddo
                  
               
                  ! get different negative integer for each ensemble member and
                  ! each domain
                  
                  do r=1,obsPert(nest,k)%N_y
                     do c=1,obsPert(nest,k)%N_x
                        domain_id = gindex(c,r)
                        call get_init_Pert_rseed(LIS_rc%nensem(nest), &
                             N_domain, LIS_rc%gnc(nest), &
                             LIS_rc%gnr(nest),  &
                             ens_id, domain_id, &
                             init_Pert_rseed (:,c,r,:))
                     enddo
                  enddo
                  deallocate(gindex)
               
               ! initialize first row of Forcepert_rseed (for first domain)
                  do n=1,LIS_rc%nensem(nest)
                     obsPert(nest,k)%obspert_rseed(1,:,:,n) = &
                          init_Pert_rseed( n, :,:,1)
                  enddo
               else
                  allocate(obsPert(nest,k)%obspert(objcount,&
                       1, 1, LIS_rc%nensem(nest)))
               endif
               ! initial call to get_pert  (for first domain)
               !
               call ESMF_StateGet(Pert_State(nest), itemCount=objcount,rc=status)
               call LIS_verify(status)
               
               allocate(pertobjs(objcount))
               allocate(pertField(objcount))
               
               call ESMF_StateGet(Pert_State(nest), itemNameList=pertobjs, rc=status)
               call LIS_verify(status)
               
               if(LIS_masterproc) then 
                  allocate(std(objcount,LIS_rc%obs_gnc(k), LIS_rc%obs_gnr(k)))
                  allocate(std1(LIS_rc%obs_ngrid(k)))
               else
                  allocate(std(1,1,1))
                  allocate(std1(LIS_rc%obs_ngrid(k)))
               endif

               do i=1,objcount
                  call ESMF_StateGet(Pert_State(nest), pertobjs(i),pertField(i),&
                       rc=status)
                  call LIS_verify(status)
                  
                  if(LIS_rc%obs_ngrid(k).gt.0) then 
                     call ESMF_AttributeGet(pertField(i),"Standard Deviation",&
                          std1,rc=status)
                     call LIS_verify(status)
                  endif

                  call LIS_gather_1dgrid_to_2dgrid_obs(nest, k,global_std, std1)
                  
                  if(LIS_masterproc) then 
                     std(i,:,:) = global_std(:,:)
                  endif

                  deallocate(global_std)

               enddo
               
               if(LIS_masterproc) then 
                  call get_pert(                                           &
                       objcount, LIS_rc%nensem(nest),&
                       LIS_rc%obs_gnc(k),LIS_rc%obs_gnr(k),  &
                       obsPert(nest,k)%N_x, obsPert(nest,k)%N_y, &
                       dx, dy, dtstep,                                     &
                       std,                              & 
                       obsPert(nest,k)%obspert_param, &
                       obsPert(nest,k)%obspert_rseed, &
                       obsPert(nest,k)%obspert_ntrmdt, &
                       obsPert(nest,k)%obspert,   &                
                       initialize_rseed=.true.,    &
                       initialize_ntrmdt=.true.            )
                  
                  deallocate(init_Pert_rseed)
                  deallocate(ens_id)
               endif
               deallocate(std)
               deallocate(std1)
            endif
         enddo
      endif

    end subroutine gmaoPert_setup
!BOP
! !ROUTINE: gmaoperturb
!
! !INTERFACE:      
    subroutine gmaoperturb(id, n, k, Base_State, Pert_State)
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
      real,             allocatable :: std(:,:,:)
      real,             allocatable :: std1(:)
      real,     allocatable         :: global_std(:,:)
      real,     allocatable         :: ltmp(:,:,:)
!      real                          :: ltmp(LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_rc%nensem(n))

      write(unit=temp,fmt='(i2.2)')n
      read(unit=temp,fmt='(2a1)') nestid

      if(id.eq.3) then 
         allocate(ltmp(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k),LIS_rc%nensem(n)))
      else
         allocate(ltmp(LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_rc%nensem(n)))
      endif

      if(.not.isXYcorrEnabled(Pert_State)) then 

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
            allocate(std(objcount,LIS_rc%obs_lnc(k), LIS_rc%obs_lnr(k)))
            allocate(std1(LIS_rc%obs_ngrid(k)))
         else
            allocate(std(objcount,LIS_rc%lnc(n), LIS_rc%lnr(n)))
            allocate(std1(LIS_rc%ngrid(n)))
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
                  call ESMF_AttributeGet(pertField(i),"Standard Deviation",std1,&
                    rc=status)
                  call LIS_verify(status)
               endif
            else
               if(LIS_rc%ngrid(n).gt.0) then 
                  call ESMF_AttributeGet(pertField(i),"Standard Deviation",std1,&
                    rc=status)
                  call LIS_verify(status)
               endif
            endif

            if(id.eq.1) then  
               std(i,:,:) = 0.0 
               do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)            
                  do m=1, LIS_rc%nensem(n)
                     tindex = (t-1)*LIS_rc%nensem(n)+m
                     
                     col = LIS_domain(n)%tile(tindex)%col
                     row = LIS_domain(n)%tile(tindex)%row
                     gid = LIS_domain(n)%gindex(col,row)
                     
                     std(i,col,row) = std(i,col,row) + std1(gid)/LIS_rc%nensem(n) 
                  enddo
               enddo
            elseif(id.eq.2) then !lsm state space
               std(i,:,:) = 0.0 
               do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)/LIS_rc%nensem(n)            
                  do m=1, LIS_rc%nensem(n)
                     tindex = (t-1)*LIS_rc%nensem(n)+m
                     
                     col = LIS_surface(n,LIS_rc%lsm_index)%tile(tindex)%col
                     row = LIS_surface(n,LIS_rc%lsm_index)%tile(tindex)%row
                     gid = LIS_domain(n)%gindex(col,row)
                     
                     std(i,col,row) = std(i,col,row) + std1(gid)/LIS_rc%nensem(n) 
                  enddo
               enddo
            elseif(id.eq.3) then 
               std(i,:,:) = 0.0 
               do t=1,LIS_rc%obs_ngrid(k)
                  col = LIS_obs_domain(n,k)%col(t)
                  row = LIS_obs_domain(n,k)%row(t)
                  std(i,col,row) = std1(t)
               enddo
            endif
         enddo
         
         
         if(id.eq.1) then 
            call LIS_getDomainResolutions(n,dx,dy)
            
            dtstep = LIS_rc%ts
            call get_pert(                                           &
                 LIS_rc%nforcepert, LIS_rc%nensem(n), &
                 LIS_rc%lnc(n), LIS_rc%lnr(n),       &
                 forcPert(n,k)%N_x, forcPert(n,k)%N_y, &
                 dx, dy, dtstep,                  &
                 std, & 
                 forcPert(n,k)%forcepert_param,  &
                 forcPert(n,k)%Forcepert_rseed,  &
                 forcPert(n,k)%Forcepert_ntrmdt, &
                 forcPert(n,k)%Forcepert)
            
            do i=1,objcount
               call ESMF_StateGet(Pert_State, pertobjs(i),pertField(i),&
                    rc=status)
               call LIS_verify(status)
               call ESMF_FieldGet(pertField(i),localDE=0,&
                    farrayPtr=pertdata1d,rc=status)
               call LIS_verify(status)
               
               do t=1,LIS_rc%ntiles(n)
                  col = LIS_domain(n)%tile(t)%col
                  row = LIS_domain(n)%tile(t)%row
                  ensem = LIS_domain(n)%tile(t)%ensem
                  pertdata1d(t) = forcpert(n,k)%forcepert(i,col,row,ensem)
               enddo
            enddo
         elseif(id.eq.2) then 

            call LIS_getDomainResolutions(n,dx,dy)
            dtstep = LIS_rc%ts
            call get_pert(                                           &
                 objcount, LIS_rc%nensem(n), LIS_rc%lnc(n), LIS_rc%lnr(n), &
                 progPert(n,k)%N_x, progPert(n,k)%N_y, &
                 dx, dy, dtstep,                  &
                 std,                           & 
                 progPert(n,k)%progpert_param,  &
                 progPert(n,k)%progpert_rseed,  &
                 progPert(n,k)%progpert_ntrmdt, &
                 progPert(n,k)%progpert)
            
            do i=1,objcount
               call ESMF_StateGet(Pert_State, pertobjs(i),pertField(i),&
                    rc=status)
               call LIS_verify(status, 'ESMF_StateGet failed in gmaoperturb')
               call ESMF_FieldGet(pertField(i),localDE=0,&
                    farrayPtr=pertdata1d,rc=status)
               call LIS_verify(status,'ESMF_FieldGet failed in gmaoperturb')
               
               do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
                  col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
                  row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
                  ensem = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%ensem
                  pertdata1d(t) = progpert(n,k)%progpert(i,col,row,ensem)
               enddo
            enddo
            
         elseif(id.eq.3) then 
            call LIS_getObsDomainResolutions(k,dx,dy)
            dtstep = LIS_rc%ts
            call get_pert(                             &
                 objcount, LIS_rc%nensem(n),           &
                 LIS_rc%obs_lnc(k), LIS_rc%obs_lnr(k), &
                 obsPert(n,k)%N_x, obsPert(n,k)%N_y,   &
                 dx, dy, dtstep,                       &
                 std,                                  & 
                 obsPert(n,k)%obspert_param,           &
                 obsPert(n,k)%obspert_rseed,           &
                 obsPert(n,k)%obspert_ntrmdt,          &
                 obsPert(n,k)%obspert)
            
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
                                obspert(n,k)%obspert(i,c,r,m)
                        enddo
                     endif
                  enddo
               enddo
            enddo
         endif
         
         deallocate(pertobjs)
         deallocate(pertField)
         deallocate(std)
         deallocate(std1)
         
      else !xy corrs enabled

         call ESMF_StateGet(Pert_State, itemCount=objcount,rc=status)
         call LIS_verify(status,'ESMF_StateGet in gmaopert_setup')
         
         call ESMF_StateGet(Base_State, itemCount=objcount_b,rc=status)
         call LIS_verify(status,'ESMF_StateGet in gmaopert_setup')
         
         allocate(pertobjs(objcount))
         allocate(pertField(objcount))
         
         allocate(baseobjs(objcount_b))
         allocate(baseField(objcount_b))
         
         call ESMF_StateGet(Base_State, itemNameList=baseobjs, rc=status)
         call LIS_verify(status,'ESMF_StateGet in gmaopert_setup')
         
         call ESMF_StateGet(Pert_State, itemNameList=pertobjs, rc=status)
         call LIS_verify(status,'ESMF_StateGet in gmaopert_setup')

         if(id.eq.3) then 
            if(LIS_masterproc) then 
               allocate(std(objcount,LIS_rc%obs_gnc(k), LIS_rc%obs_gnr(k)))
               allocate(std1(LIS_rc%obs_ngrid(k)))
            else
               allocate(std(1,1,1))
               allocate(std1(LIS_rc%obs_ngrid(k)))
            endif
         else
            if(LIS_masterproc) then 
               allocate(std(objcount,LIS_rc%gnc(n), LIS_rc%gnr(n)))
               allocate(std1(LIS_rc%ngrid(n)))
            else
               allocate(std(1,1,1))
               allocate(std1(LIS_rc%ngrid(n)))
            endif
         endif
         do i=1,objcount
! The forcing state and perturbations state have different members. For
! lsm and obs objects, the number of elements is the same in both 
! the LSM_state (Obs_state) and LSM_pert_state (Obs_pert_state). So 
! here we need to handle these separately. 

            if(id.eq.1) then 
               baseobjs(i) = pertobjs(i)
            endif
            call ESMF_StateGet(Base_State, trim(baseobjs(i)),baseField(i),&
                 rc=status)
            call LIS_verify(status,'ESMF_StateGet in gmaopert_setup')
            call ESMF_FieldGet(baseField(i),localDE=0,farrayPtr=basedata,rc=status)
            call LIS_verify(status,'ESMF_FieldGet in gmaopert_setup')
            
            call ESMF_StateGet(Pert_State, trim(pertobjs(i)),pertField(i),&
                 rc=status)
            call LIS_verify(status,'ESMF_StateGet in gmaopert_setup')
            
            if(id.eq.3) then 
               if(LIS_rc%obs_ngrid(k).gt.0) then 
                  call ESMF_AttributeGet(pertField(i),"Standard Deviation",std1,&
                       rc=status)
                  call LIS_verify(status,'ESMF_AttributeGet (1) failed in gmaopert_setup')
               endif
               call LIS_gather_1dgrid_to_2dgrid_obs(n, k, global_std, std1)
               if(LIS_masterproc) then 
                  std(i,:,:) = global_std(:,:) 
               endif
               deallocate(global_std)
            else         
               if(LIS_rc%ngrid(n).gt.0) then 
                  call ESMF_AttributeGet(pertField(i),"Standard Deviation",std1,&
                       rc=status)
                  call LIS_verify(status,'ESMF_AttributeGet (2)failed in gmaopert_setup')
               endif
               
               call LIS_gather_1dgrid_to_2dgrid(n, global_std, std1)

               if(LIS_masterproc) then 
                  std(i,:,:) = global_std(:,:) 
               endif
               deallocate(global_std)
            endif

         enddo
         
         if(id.eq.1) then 
            if(LIS_masterproc) then 
               call LIS_getDomainResolutions(n,dx,dy)
               
               dtstep = LIS_rc%ts
               call get_pert(                                           &
                    LIS_rc%nforcepert, LIS_rc%nensem(n), &
                    LIS_rc%gnc(n), LIS_rc%gnr(n),       &
                    forcPert(n,k)%N_x, forcPert(n,k)%N_y, &
                    dx, dy, dtstep,                  &
                    std, & 
                    forcPert(n,k)%forcepert_param,  &
                    forcPert(n,k)%Forcepert_rseed,  &
                    forcPert(n,k)%Forcepert_ntrmdt, &
                    forcPert(n,k)%Forcepert)
            endif
            do i=1,objcount
               do m=1,LIS_rc%nensem(n)
                  call LIS_scatter_global_to_local_grid(n, &
                       forcPert(n,k)%ForcePert(i,:,:,m),ltmp(:,:,m))
               enddo

               call ESMF_StateGet(Pert_State, pertobjs(i),pertField(i),&
                    rc=status)
               call LIS_verify(status)
               call ESMF_FieldGet(pertField(i),localDE=0,&
                    farrayPtr=pertdata1d,rc=status)
               call LIS_verify(status)
               
               do t=1,LIS_rc%ntiles(n)
                  col = LIS_domain(n)%tile(t)%col
                  row = LIS_domain(n)%tile(t)%row
                  ensem = LIS_domain(n)%tile(t)%ensem
                  pertdata1d(t) = ltmp(col,row,ensem)
               enddo
            enddo
         elseif(id.eq.2) then 
            if(LIS_masterproc) then 
               call LIS_getDomainResolutions(n,dx,dy)
               dtstep = LIS_rc%ts
               call get_pert(                                           &
                    objcount, LIS_rc%nensem(n), LIS_rc%gnc(n), LIS_rc%gnr(n), &
                    progPert(n,k)%N_x, progPert(n,k)%N_y, &
                    dx, dy, dtstep,                  &
                    std,                           & 
                    progPert(n,k)%progpert_param,  &
                    progPert(n,k)%progpert_rseed,  &
                    progPert(n,k)%progpert_ntrmdt, &
                    progPert(n,k)%progpert)
            endif
            do i=1,objcount
               do m=1,LIS_rc%nensem(n)

                  call LIS_scatter_global_to_local_grid(n, &
                       progPert(n,k)%progPert(i,:,:,m),ltmp(:,:,m))
               enddo
               call ESMF_StateGet(Pert_State, pertobjs(i),pertField(i),&
                    rc=status)
               call LIS_verify(status)
               call ESMF_FieldGet(pertField(i),localDE=0,&
                    farrayPtr=pertdata1d,rc=status)
               call LIS_verify(status)
               
               do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
                  col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
                  row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
                  ensem = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%ensem
                  pertdata1d(t) = ltmp(col,row,ensem)
               enddo
            enddo
            
         elseif(id.eq.3) then 
            if(LIS_masterproc) then 
               call LIS_getObsDomainResolutions(k,dx,dy)
               dtstep = LIS_rc%ts
               call get_pert(                                           &
                    objcount, LIS_rc%nensem(n), &
                    LIS_rc%obs_gnc(k), LIS_rc%obs_gnr(k),&
                    obsPert(n,k)%N_x, obsPert(n,k)%N_y, &
                    dx, dy, dtstep,                  &
                    std,                           & 
                    obsPert(n,k)%obspert_param,  &
                    obsPert(n,k)%obspert_rseed,  &
                    obsPert(n,k)%obspert_ntrmdt, &
                    obsPert(n,k)%obspert)
            endif

            do i=1,objcount            
               do m =1, LIS_rc%nensem(n)
                  call LIS_scatter_global_to_local_grid_obs(n,k,&
                       obsPert(n,k)%obsPert(i,:,:,m),ltmp(:,:,m))
               enddo
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
                                ltmp(c,r,m)
                        enddo
                     endif
                  enddo
               enddo
            enddo
         endif
         
         deallocate(pertobjs)
         deallocate(pertField)
         deallocate(std)
         deallocate(std1)
         deallocate(ltmp)

      endif
    end subroutine gmaoperturb


!BOP
! 
! !ROUTINE: gmaopert_writerestart
!  \label{gmaopert_writerestart}
! 
! !INTERFACE: 
    subroutine gmaopert_writerestart(n)
! !USES: 

! !ARGUMENTS: 
      implicit none

      integer, intent(in)      :: n 
! !DESCRIPTION: 
! 
! Write one single file that includes the restart info for forcing, states 
! and observations

!EOP      
      integer                  :: k, kk,i, col, row, ensem, t
      real, allocatable        :: pertdata1d(:)
      real, allocatable        :: pertdata1d_patch(:)
      integer, allocatable     :: pertdata1d_int(:)
      integer, allocatable     :: pertdata1d_patch_int(:)
      character*100            :: filen
      integer                  :: ftn 

      if ( LIS_masterproc ) then
         
         call LIS_create_output_directory('DAPERT')         
         call LIS_create_dapert_filename(n,filen)
         
         ftn = LIS_getNextUnitNumber()
         write(LIS_logunit,*) '[INFO] Writing Perturbations restart ', &
              trim(filen)
         open(ftn, file = trim(filen), form='unformatted')
         
      endif
      if(LIS_rc%perturb_forcing .ne."none") then 
         if(.not. f_xyCorr) then 
            allocate(pertdata1d(LIS_rc%ntiles(n)))
            
            k = 1
            do i=1,LIS_rc%nforcepert
               
               do t=1,LIS_rc%ntiles(n)
                  col = LIS_domain(n)%tile(t)%col
                  row = LIS_domain(n)%tile(t)%row
                  ensem = LIS_domain(n)%tile(t)%ensem
                  pertdata1d(t) = forcpert(n,k)%forcepert_ntrmdt(i,col,row,ensem)
               enddo
               
               call LIS_writevar_restart(ftn,n,pertdata1d)         
            enddo
            deallocate(pertdata1d)
            allocate(pertdata1d_int(LIS_rc%ntiles(n)))
            
            do i=1,NRANDSEED            
               do t=1,LIS_rc%ntiles(n)
                  col = LIS_domain(n)%tile(t)%col
                  row = LIS_domain(n)%tile(t)%row
                  ensem = LIS_domain(n)%tile(t)%ensem
                  pertdata1d_int(t) = &
                       forcpert(n,k)%forcepert_rseed(i,col,row,ensem)
               enddo
               call LIS_writevar_restart(ftn,n,pertdata1d_int)         
            enddo
            deallocate(pertdata1d_int)

         else
            if(LIS_masterproc) then 
               k = 1
               do i=1,LIS_rc%nforcepert               
                  write(ftn) forcpert(n,k)%forcepert_ntrmdt(i,:,:,:)
               enddo
               
               do i=1,NRANDSEED
                  write(ftn) forcpert(n,k)%forcepert_rseed(i,:,:,:)
               enddo
            endif
         endif
      endif

      if(.not. p_xyCorr) then 

         do k=1,LIS_rc%nperts
            allocate(pertdata1d_patch(LIS_rc%npatch(n,LIS_rc%lsm_index)))

            if(trim(LIS_rc%perturb_state(k)).ne."none") then 
               do i =1, LIS_rc%nstvars(k)
                  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
                     col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
                     row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
                     ensem = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%ensem
                     pertdata1d_patch(t) =  &
                          progpert(n,k)%progpert_ntrmdt(i,col,row,ensem)
                  enddo
                  
                  call LIS_writevar_restart(ftn,n,1,pertdata1d_patch)   
               enddo
            endif

            deallocate(pertdata1d_patch)
            allocate(pertdata1d_patch_int(LIS_rc%npatch(n,LIS_rc%lsm_index)))

            if(trim(LIS_rc%perturb_state(k)).ne."none") then 
               do i =1, NRANDSEED
                  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
                     col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
                     row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
                     ensem = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%ensem
                     pertdata1d_patch_int(t) =  &
                          progpert(n,k)%progpert_rseed(i,col,row,ensem)
                  enddo
                  
                  call LIS_writevar_restart(ftn,n,1,pertdata1d_patch_int)   
               enddo
            endif

            deallocate(pertdata1d_patch_int)
         enddo

      else
         if(LIS_masterproc) then 
            do k=1,LIS_rc%nperts
               if(trim(LIS_rc%perturb_state(k)).ne."none") then 
                  do i =1, LIS_rc%nstvars(k)
                     write(ftn) &
                          progpert(n,k)%progpert_ntrmdt(i,:,:,:)
                  enddo

                  do i =1, NRANDSEED
                     write(ftn) &
                          progpert(n,k)%progpert_rseed(i,:,:,:)
                  enddo
               endif
            enddo
         endif
      endif

      if(.not. o_xyCorr) then 
         do k=1,LIS_rc%ndas
            if(trim(LIS_rc%perturb_obs(k)).ne."none") then 
               allocate(pertdata1d(LIS_rc%obs_ngrid(k)*LIS_rc%nensem(n)))
               do i=1, nobs(n,k)
                  do t=1,LIS_rc%obs_ngrid(k)
                     col = LIS_obs_domain(n,k)%col(t)
                     row = LIS_obs_domain(n,k)%row(t)
                     do kk=1,LIS_rc%nensem(n)
                        pertdata1d(kk+(t-1)*LIS_rc%nensem(n)) =  &
                             obspert(n,k)%obspert_ntrmdt(i,col,row,kk)
                     enddo
                  enddo

                  call writevar_obspert_restart(ftn,n,k,pertdata1d)   
               enddo
               deallocate(pertdata1d)

               allocate(pertdata1d_int(LIS_rc%obs_ngrid(k)*LIS_rc%nensem(n)))
               do i=1, NRANDSEED
                  do t=1,LIS_rc%obs_ngrid(k)
                     col = LIS_obs_domain(n,k)%col(t)
                     row = LIS_obs_domain(n,k)%row(t)
                     do kk=1,LIS_rc%nensem(n)
                        pertdata1d_int(kk+(t-1)*LIS_rc%nensem(n)) =  &
                             obspert(n,k)%obspert_rseed(i,col,row,kk)
                     enddo
                  enddo

                  call writevar_obspert_restart_int(ftn,n,k,pertdata1d_int)   
               enddo
               deallocate(pertdata1d_int)

            endif
         enddo
      else
         if(LIS_masterproc) then 
            do k=1,LIS_rc%ndas
               if(trim(LIS_rc%perturb_obs(k)).ne."none") then 
                  do i=1, nobs(n,k)
                     write(ftn) obspert(n,k)%obspert_ntrmdt(i,:,:,:)
                  enddo
                  do i=1, NRANDSEED
                     write(ftn) obspert(n,k)%obspert_rseed(i,:,:,:)
                  enddo
               endif
            enddo
         endif
      endif
      if(LIS_masterproc) then 
         call LIS_releaseUnitNumber(ftn)
         write(LIS_logunit,*) '[INFO] Done writing Perturbations restart ', &
              trim(filen)
      endif

    end subroutine gmaopert_writerestart

!BOP
! 
! !ROUTINE: gmaopert_readrestart
! \label{gmaopert_readrestart}
! 
! !INTERFACE: 
    subroutine gmaopert_readrestart()
! !USES: 

      implicit none
! !ARGUMENTS: 

! !DESCRIPTION: 
!   This routine reads from a perturbations restart file
!EOP
      integer                   :: n 
      integer                   :: k, kk,i,t, ftn, col, row, ensem
      integer                   :: yr,mo,da,hr,mn,ss, doy
      integer                   :: status
      character*100             :: filen
      real*8                    :: time
      real                      :: gmt
      real, allocatable         :: pertdata1d(:)
      real, allocatable         :: pertdata1d_obs(:)
      real, allocatable         :: pertdata1d_patch(:)
      integer, allocatable      :: pertdata1d_int(:)
      integer, allocatable      :: pertdata1d_obs_int(:)
      integer, allocatable      :: pertdata1d_patch_int(:)
      real, allocatable         :: dummy_var(:,:,:)

      do n = 1, LIS_rc%nnest

         ftn = LIS_getNextUnitNumber()

         if(LIS_rc%runmode.eq."ensemble smoother") then 
            if(LIS_rc%iterationid(n).gt.1) then 
               if(LIS_rc%pertrestartInterval.eq.2592000) then
                  call ESMF_TimeGet(LIS_twStartTime,yy=yr,mm=mo,&
                       dd=da,calendar=LIS_calendar,rc=status)
                  hr = 0 
                  mn = 0 
                  ss = 0 
                  call LIS_tick(time,doy,gmt,yr,mo,da,hr,mn,ss,(-1)*LIS_rc%ts)
               else
                  call ESMF_TimeGet(LIS_twStartTime,yy=yr,mm=mo,&
                       dd=da,calendar=LIS_calendar,rc=status)
                  hr = 0 
                  mn = 0 
                  ss = 0 
               endif
               call LIS_create_dapert_filename(n,filen, yr, mo, da, hr, mn, ss)
               LIS_rc%pertRestartFile(n) = filen
            endif
         endif
         
         open(ftn,file=trim(LIS_rc%pertRestartFile(n)), form='unformatted')
         write(LIS_logunit,*) '[INFO] Reading perturbations restart file ',&
              trim(LIS_rc%pertRestartFile(n))
         
         if(LIS_rc%perturb_forcing .ne."none") then 
            if(.not. f_xyCorr) then 
               allocate(pertdata1d(LIS_rc%ntiles(n)))
               
               k = 1
               do i=1,LIS_rc%nforcepert
                  
                  call LIS_readvar_restart(ftn,n,pertdata1d)
                  
                  do t=1,LIS_rc%ntiles(n)
                     col = LIS_domain(n)%tile(t)%col
                     row = LIS_domain(n)%tile(t)%row
                     ensem = LIS_domain(n)%tile(t)%ensem
                     forcpert(n,k)%forcepert_ntrmdt(i,col,row,ensem) = pertdata1d(t)
                  enddo
               enddo
               
               deallocate(pertdata1d)
               
               allocate(pertdata1d_int(LIS_rc%ntiles(n)))
               
               do i=1,NRANDSEED
                  
                  call LIS_readvar_restart(ftn,n,pertdata1d_int)
                  
                  do t=1,LIS_rc%ntiles(n)
                     col = LIS_domain(n)%tile(t)%col
                     row = LIS_domain(n)%tile(t)%row
                     ensem = LIS_domain(n)%tile(t)%ensem
                     forcpert(n,k)%forcepert_rseed(i,col,row,ensem) = pertdata1d_int(t)
                  enddo
               enddo
               deallocate(pertdata1d_int)
               
            else
               if(LIS_masterproc) then
                  k = 1
                  do i=1,LIS_rc%nforcepert                  
                     read(ftn) forcpert(n,k)%forcepert_ntrmdt(i,:,:,:)
                  enddo
                  do i=1,NRANDSEED
                     read(ftn) forcpert(n,k)%forcepert_rseed(i,:,:,:)
                  enddo
               else
! to keep the binary reading in sync when there is a mix of file reading that
! happens only on one processor and across processors. 
                  allocate(dummy_var(&
                       LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%nensem(n)))
                  do i =1, LIS_rc%nforcepert
                     read(ftn) dummy_var
                  enddo
                  do i =1, NRANDSEED
                     read(ftn) dummy_var
                  enddo
                  deallocate(dummy_var)
               endif
               
            endif
         endif
         if(.not. p_xyCorr) then 
            do k=1, LIS_rc%nperts

               if(trim(LIS_rc%perturb_state(k)).ne."none") then 
                  allocate(pertdata1d_patch(LIS_rc%npatch(n,LIS_rc%lsm_index)))
                  
                  do i =1, LIS_rc%nstvars(k)                  
                     call LIS_readvar_restart(ftn,n,1,pertdata1d_patch)   
                     
                     do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
                        col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
                        row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
                        ensem = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%ensem
                        progpert(n,k)%progpert_ntrmdt(i,col,row,ensem) = &
                             pertdata1d_patch(t) 
                     enddo
                     
                  enddo
                  deallocate(pertdata1d_patch)
                  
                  allocate(pertdata1d_patch_int(LIS_rc%npatch(n,LIS_rc%lsm_index)))

                  do i =1, NRANDSEED
                     call LIS_readvar_restart(ftn,n,1,pertdata1d_patch_int)   
                     
                     do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
                        col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
                        row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
                        ensem = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%ensem
                        progpert(n,k)%progpert_rseed(i,col,row,ensem) = &
                             pertdata1d_patch_int(t) 
                     enddo
                     
                  enddo
                  deallocate(pertdata1d_patch_int)

               endif
            enddo
         else
            if(LIS_masterproc) then 
               do k=1, LIS_rc%nperts
                  if(trim(LIS_rc%perturb_state(k)).ne."none") then 
                  
                     do i =1, LIS_rc%nstvars(k)                  
                        read(ftn) progPert(n,k)%progpert_ntrmdt(i,:,:,:)
                     enddo
                     do i =1, NRANDSEED
                        read(ftn) progPert(n,k)%progpert_rseed(i,:,:,:)
                     enddo
                  endif
               enddo
            else
! to keep the binary reading in sync when there is a mix of file reading that
! happens only on one processor and across processors. 
               allocate(dummy_var(&
                    LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_rc%nensem(n)))
               do k=1, LIS_rc%nperts
                  if(trim(LIS_rc%perturb_state(k)).ne."none") then 
                     
                     do i =1, LIS_rc%nstvars(k)                  
                        read(ftn) dummy_var
                     enddo
                     do i =1, NRANDSEED
                        read(ftn) dummy_var
                     enddo
                  endif
               enddo
               deallocate(dummy_var)
            endif
         endif   

         if(.not. o_xyCorr) then       
            do k=1, LIS_rc%ndas

               if(trim(LIS_rc%perturb_obs(k)).ne."none") then 
                  allocate(pertdata1d_obs(LIS_rc%obs_ngrid(k)*LIS_rc%nensem(n)))

                  do i=1, nobs(n,k)
                     call readvar_obspert_restart(ftn,n,k,pertdata1d_obs)   

                     do t=1,LIS_rc%obs_ngrid(k)
                        col = LIS_obs_domain(n,k)%col(t)
                        row = LIS_obs_domain(n,k)%row(t)
                        do kk=1,LIS_rc%nensem(n)
                           obspert(n,k)%obspert_ntrmdt(i,col,row,kk) = & 
                                pertdata1d_obs(kk+(t-1)*LIS_rc%nensem(n))
                               
                        enddo
                     enddo
                  enddo
                  deallocate(pertdata1d_obs)
                  
                  allocate(pertdata1d_obs_int(LIS_rc%obs_ngrid(k)*LIS_rc%nensem(n)))
                  

                  do i=1, NRANDSEED
                     call readvar_obspert_restart_int(ftn,n,k,pertdata1d_obs_int)   

                     do t=1,LIS_rc%obs_ngrid(k)
                        col = LIS_obs_domain(n,k)%col(t)
                        row = LIS_obs_domain(n,k)%row(t)
                        do kk=1,LIS_rc%nensem(n)
                           obspert(n,k)%obspert_rseed(i,col,row,kk) = & 
                                pertdata1d_obs_int(kk+(t-1)*LIS_rc%nensem(n))
                               
                        enddo
                     enddo
                  enddo
               endif
               deallocate(pertdata1d_obs_int)

            enddo
         else
            if(LIS_masterproc) then 
               do k=1, LIS_rc%ndas
                  if(trim(LIS_rc%perturb_obs(k)).ne."none") then 
                     do i=1, nobs(n,k)
                        read(ftn) obsPert(n,k)%obspert_ntrmdt(i,:,:,:)
                     enddo
                     do i=1, NRANDSEED
                        read(ftn) obsPert(n,k)%obspert_rseed(i,:,:,:)
                     enddo
                  end if
               enddo
            else
               do k=1, LIS_rc%nperts
                  allocate(dummy_var(&
                       LIS_rc%obs_gnc(k),LIS_rc%obs_gnr(k),LIS_rc%nensem(n)))
                  if(trim(LIS_rc%perturb_obs(k)).ne."none") then  
                     do i =1, nobs(n,k)
                        read(ftn) dummy_var
                     enddo
                     do i =1, NRANDSEED
                        read(ftn) dummy_var
                     enddo
                  endif
                  deallocate(dummy_var)

               enddo

            endif
         endif
         
         call LIS_releaseUnitNumber(ftn)
         write(LIS_logunit,*) '[INFO] Finished reading perturbations restart file ',&
              trim(LIS_rc%pertRestartFile(n))

      enddo
    end subroutine gmaopert_readrestart

!BOP
! !ROUTINE: writevar_obspert_restart
! \label{writevar_obspert_restart}
!
! !INTERFACE:
  subroutine writevar_obspert_restart(ftn, n, k, var)
! !USES:
!    use LIS_DAobservationsMod

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: ftn
    integer, intent(in) :: n
    integer, intent(in) :: k
    real, intent(in)    :: var(LIS_rc%obs_ngrid(k)*LIS_rc%nensem(n))
!
! !DESCRIPTION:
!  Writes the observation perturbation to a binary file in a gridded format. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the observation space
!  \end{description}
!EOP

    real, allocatable :: gtmp(:)
    real, allocatable :: gtmp1(:)
    integer           :: ierr
    integer           :: odeltas
    integer           :: stid,tid
    integer           :: i,c,r,l,t,ntiles,gid,count1

    if(LIS_masterproc) then 
       allocate(gtmp1(LIS_rc%obs_glbngrid(k)*LIS_rc%nensem(n)))
    else
       allocate(gtmp1(1))
    endif
#if (defined SPMD)
    odeltas =  LIS_obsens_gdeltas(k,LIS_localPet)
    call MPI_GATHERV(var,odeltas,&
         MPI_REAL,gtmp1,LIS_obsens_gdeltas(k,:),&
         LIS_obsens_goffsets(k,:),&
         MPI_REAL,0,LIS_mpi_comm,ierr)
#else 
!    gtmp = var
    gtmp1 = var
#endif
    if(LIS_masterproc) then 

#if 0 
       do l=1,LIS_npes
          do r=LIS_nss_obs_halo_ind(k,l),LIS_nse_obs_halo_ind(k,l)
             do c=LIS_ews_obs_halo_ind(k,l),LIS_ewe_obs_halo_ind(k,l)
                gid = c+(r-1)*LIS_rc%obs_gnc(k)
                ntiles = LIS_obs_domain(n,k)%ntiles_pergrid(gid)
                stid = LIS_obs_domain(n,k)%str_tind(gid)
                if(ntiles.ne.0) then
                   if(r.ge.LIS_nss_obs_ind(k,l).and.&
                        r.le.LIS_nse_obs_ind(k,l).and.&
                        c.ge.LIS_ews_obs_ind(k,l).and.&
                        c.le.LIS_ewe_obs_ind(k,l)) then !points not in halo  
                      do t=1,ntiles*LIS_rc%nensem(n)
                         tid = (stid - 1)*LIS_rc%nensem(n) + t
                         gtmp(tid) = gtmp1(count1)
                         count1 = count1 + 1
                      enddo
                   endif
!                   count1 = count1 + LIS_rc%nensem(n)
                endif
             enddo
          enddo
       enddo
#endif
       write(ftn) gtmp1
!       deallocate(gtmp)
    endif
    deallocate(gtmp1)
  end subroutine Writevar_obspert_restart

!BOP
! !ROUTINE: writevar_obspert_restart_int
! \label{writevar_obspert_restart_int}
!
! !INTERFACE:
  subroutine writevar_obspert_restart_int(ftn, n, k, var)
! !USES:
!    use LIS_DAobservationsMod

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: ftn
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in)    :: var(LIS_rc%obs_ngrid(k)*LIS_rc%nensem(n))
!
! !DESCRIPTION:
!  Writes the observation perturbation to a binary file in a gridded format. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the observation space
!  \end{description}
!EOP

    integer, allocatable :: gtmp(:)
    integer, allocatable :: gtmp1(:)
    integer           :: ierr
    integer           :: odeltas
    integer           :: stid,tid
    integer           :: i,c,r,l,t,ntiles,gid,count1

    if(LIS_masterproc) then 
       allocate(gtmp1(LIS_rc%obs_glbngrid(k)*LIS_rc%nensem(n)))
    else
       allocate(gtmp1(1))
    endif
#if (defined SPMD)
    odeltas =  LIS_obsens_gdeltas(k,LIS_localPet)
    call MPI_GATHERV(var,odeltas,&
         MPI_INTEGER,gtmp1,LIS_obsens_gdeltas(k,:),&
         LIS_obsens_goffsets(k,:),&
         MPI_INTEGER,0,LIS_mpi_comm,ierr)
#else 
!    gtmp = var
    gtmp1 = var
#endif
    if(LIS_masterproc) then 

#if 0 
       do l=1,LIS_npes
          do r=LIS_nss_obs_halo_ind(k,l),LIS_nse_obs_halo_ind(k,l)
             do c=LIS_ews_obs_halo_ind(k,l),LIS_ewe_obs_halo_ind(k,l)
                gid = c+(r-1)*LIS_rc%obs_gnc(k)
                ntiles = LIS_obs_domain(n,k)%ntiles_pergrid(gid)
                stid = LIS_obs_domain(n,k)%str_tind(gid)
                if(ntiles.ne.0) then
                   if(r.ge.LIS_nss_obs_ind(k,l).and.&
                        r.le.LIS_nse_obs_ind(k,l).and.&
                        c.ge.LIS_ews_obs_ind(k,l).and.&
                        c.le.LIS_ewe_obs_ind(k,l)) then !points not in halo  
                      do t=1,ntiles*LIS_rc%nensem(n)
                         tid = (stid - 1)*LIS_rc%nensem(n) + t
                         gtmp(tid) = gtmp1(count1)
                         count1 = count1 + 1
                      enddo
                   endif
!                   count1 = count1 + LIS_rc%nensem(n)
                endif
             enddo
          enddo
       enddo
#endif
       write(ftn) gtmp1
!       deallocate(gtmp)
    endif
    deallocate(gtmp1)
  end subroutine Writevar_obspert_restart_int


!BOP
! !ROUTINE: readvar_obspert_restart
! \label{readvar_obspert_restart}
!
! !INTERFACE:
  subroutine readvar_obspert_restart(ftn, n, k, var)
! !USES:
    use LIS_DAobservationsMod

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: ftn
    integer, intent(in) :: n
    integer, intent(in) :: k
    real                :: var(LIS_rc%obs_ngrid(k)*LIS_rc%nensem(n))
!
! !DESCRIPTION:
!  Writes the observation perturbation to a binary file in a gridded format. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the observation space
!  \end{description}
!EOP

    real, allocatable :: gtmp(:)
    integer           :: count1
    integer           :: c,r,t,gid,tid,ntiles,stid,enid

    allocate(gtmp(LIS_rc%obs_glbngrid(k)*LIS_rc%nensem(n)))
    read(ftn) gtmp

    stid = LIS_obs_goffsets(k,LIS_localPet)*LIS_rc%nensem(n) + 1
    enid = stid+LIS_rc%obs_ngrid(k)*LIS_rc%nensem(n) -1

    var(:) = gtmp(stid:enid)

#if 0 
    count1=1
    do r=LIS_nss_obs_halo_ind(k,LIS_localPet+1),&
         LIS_nse_obs_halo_ind(k,LIS_localPet+1)
       do c=LIS_ews_obs_halo_ind(k,LIS_localPet+1),&
            LIS_ewe_obs_halo_ind(k,LIS_localPet+1)
          gid = c+(r-1)*LIS_rc%obs_gnc(k)
          ntiles = LIS_obs_domain(n,k)%ntiles_pergrid(gid)
          stid = LIS_obs_domain(n,k)%str_tind(gid)
          do t=1,ntiles*LIS_rc%nensem(n)
             tid = (stid-1)*LIS_rc%nensem(n) + t
             var(count1) = gtmp(tid) 
             count1 = count1 + 1
          enddo
       enddo
    enddo
#endif
    deallocate(gtmp) 
  end subroutine Readvar_obspert_restart


!BOP
! !ROUTINE: readvar_obspert_restart_int
! \label{readvar_obspert_restart_int}
!
! !INTERFACE:
  subroutine readvar_obspert_restart_int(ftn, n, k, var)
! !USES:
    use LIS_DAobservationsMod

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: ftn
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer                :: var(LIS_rc%obs_ngrid(k)*LIS_rc%nensem(n))
!
! !DESCRIPTION:
!  Writes the observation perturbation to a binary file in a gridded format. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the observation space
!  \end{description}
!EOP

    integer, allocatable :: gtmp(:)
    integer           :: count1
    integer           :: c,r,t,gid,tid,ntiles,stid,enid

    allocate(gtmp(LIS_rc%obs_glbngrid(k)*LIS_rc%nensem(n)))
    read(ftn) gtmp

    stid = LIS_obs_goffsets(k,LIS_localPet)*LIS_rc%nensem(n) + 1
    enid = stid+LIS_rc%obs_ngrid(k)*LIS_rc%nensem(n) -1

    var(:) = gtmp(stid:enid)

#if 0 
    count1=1
    do r=LIS_nss_obs_halo_ind(k,LIS_localPet+1),&
         LIS_nse_obs_halo_ind(k,LIS_localPet+1)
       do c=LIS_ews_obs_halo_ind(k,LIS_localPet+1),&
            LIS_ewe_obs_halo_ind(k,LIS_localPet+1)
          gid = c+(r-1)*LIS_rc%obs_gnc(k)
          ntiles = LIS_obs_domain(n,k)%ntiles_pergrid(gid)
          stid = LIS_obs_domain(n,k)%str_tind(gid)
          do t=1,ntiles*LIS_rc%nensem(n)
             tid = (stid-1)*LIS_rc%nensem(n) + t
             var(count1) = gtmp(tid) 
             count1 = count1 + 1
          enddo
       enddo
    enddo
#endif
    deallocate(gtmp) 
  end subroutine Readvar_obspert_restart_int


end module gmaopert_Mod

