!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !MODULE: S1_sigmaVVSM_Mod
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle S1 backscatter retrievals. 
!   
! !REVISION HISTORY: 
!  29 Aug 2019   Hans Lievens;   Initial Specification for sigma depth
!  9 Mar 2021    Isis Brangers, Michel Bechtold;  Adaptation for backscatter
! 29 Jun 2022:  Louise Busschaert; getting domain dims from files
! 
module S1_sigmaVVSM_Mod
! !USES: 
  use ESMF
!EOP
  implicit none
  
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: S1_sigmaVVSM_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: S1_sigma_struc

  type, public ::  S1_sigma_dec
     logical             :: startMode
     real                :: ssdev
     integer             :: sigmafield
     integer             :: nc,nr
     real,    allocatable    :: sigma(:,:)
     real,    allocatable    :: sigmatime(:,:)
  end type S1_sigma_dec

  type(S1_sigma_dec), allocatable :: S1_sigma_struc(:)

contains
!BOP
! 
! !ROUTINE: S1_sigmaVVSM_setup
! \label{S1_sigmaVVSM_setup}
! 
! !INTERFACE: 
  subroutine S1_sigmaVVSM_setup(k, OBS_State, OBS_Pert_State)
! !USES: 

    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_historyMod
    use LIS_perturbMod
    use LIS_DAobservationsMod
    use LIS_logMod
    use netcdf
   
    implicit none 

! !ARGUMENTS: 
    integer                ::  k 
    type(ESMF_State)       ::  OBS_State(LIS_rc%nnest)
    type(ESMF_State)       ::  OBS_Pert_State(LIS_rc%nnest)
! 
! !DESCRIPTION: 
!   
!   This routine completes the runtime initializations and 
!   creation of data structures required for S1 backscatter data. 
!
!   The arguments are: 
!   \begin{description}
!    \item[OBS\_State]   observation state 
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP
    integer                ::  n
    integer                ::  ftn
    integer                ::  i
    integer                ::  status
    type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  pertArrSpec
    character*100          ::  S1sigmaobsdir
    character*80           ::  S1_firstfile
    character*80           ::  S1_filename
    character*100          ::  temp
    real,  allocatable         ::  obsstd(:)
    character*1            ::  vid(2)
    character*40, allocatable  ::  vname(:)
    real,         allocatable  ::  varmin(:)
    real,         allocatable  ::  varmax(:)
    type(pert_dec_type)    :: obs_pert
    real, pointer          :: obs_temp(:,:)
    real, allocatable          :: ssdev(:)
    real                   :: dx, dy
    integer                :: NX, NY 
    integer                :: ncid
    integer                :: flist
    integer                :: ios
    character*100          :: infile 
    character*100          :: xname, yname 


    allocate(S1_sigma_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"S1 backscatter data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,S1sigmaobsdir,&
            rc=status)
       call LIS_verify(status,'S1 backscatter data directory: not defined')
       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            S1sigmaobsdir, rc=status)
       call LIS_verify(status)
    enddo

    do n=1,LIS_rc%nnest
       call ESMF_AttributeSet(OBS_State(n),"Data Update Status",&
            .false., rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(OBS_State(n),"Data Update Time",&
            -99.0, rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(OBS_State(n),"Data Assimilate Status",&
            .false., rc=status)
       call LIS_verify(status)
       
       call ESMF_AttributeSet(OBS_State(n),"Number Of Observations",&
            LIS_rc%obs_ngrid(k),rc=status)
       call LIS_verify(status)
       
    enddo

    write(LIS_logunit,*)'[INFO] read S1 backscatter data specifications'

!----------------------------------------------------------------------------
!   Create the array containers that will contain the observations and
!   the perturbations. Since there is only one layer
!   being assimilated, the array size is LIS_rc%ngrid(n). 
!----------------------------------------------------------------------------

    do n=1,LIS_rc%nnest
       
       write(unit=temp,fmt='(i2.2)') 1
       read(unit=temp,fmt='(2a1)') vid

       obsField(n) = ESMF_FieldCreate(grid=LIS_obsvecGrid(n,k),&
            arrayspec=realarrspec,&
            name="Observation"//vid(1)//vid(2), rc=status)
       call LIS_verify(status)
       
!Perturbations State
       write(LIS_logunit,*) '[INFO] Opening attributes for observations ',&
            trim(LIS_rc%obsattribfile(k))
       ftn = LIS_getNextUnitNumber()
       open(ftn,file=trim(LIS_rc%obsattribfile(k)),status='old')
       read(ftn,*)
       read(ftn,*) LIS_rc%nobtypes(k)
       read(ftn,*)
    
       allocate(vname(LIS_rc%nobtypes(k)))
       allocate(varmax(LIS_rc%nobtypes(k)))
       allocate(varmin(LIS_rc%nobtypes(k)))
       
       do i=1,LIS_rc%nobtypes(k)
          read(ftn,fmt='(a40)') vname(i)
          read(ftn,*) varmin(i),varmax(i)
          write(LIS_logunit,*) '[INFO] ',vname(i),varmin(i),varmax(i)
       enddo
       call LIS_releaseUnitNumber(ftn)  
       
       call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
       call LIS_verify(status)

       allocate(ssdev(LIS_rc%obs_ngrid(k)))

       if(trim(LIS_rc%perturb_obs(k)).ne."none") then 

          allocate(obs_pert%vname(1))
          allocate(obs_pert%perttype(1))
          allocate(obs_pert%ssdev(1))
          allocate(obs_pert%stdmax(1))
          allocate(obs_pert%zeromean(1))
          allocate(obs_pert%tcorr(1))
          allocate(obs_pert%xcorr(1))
          allocate(obs_pert%ycorr(1))
          allocate(obs_pert%ccorr(1,1))

          call LIS_readPertAttributes(1,LIS_rc%obspertAttribfile(k),&
               obs_pert)

          ssdev = obs_pert%ssdev(1)
          S1_sigma_struc(n)%ssdev =obs_pert%ssdev(1) 

          pertField(n) = ESMF_FieldCreate(arrayspec=pertArrSpec,&
               grid=LIS_obsEnsOnGrid(n,k),name="Observation"//vid(1)//vid(2),&
               rc=status)
          call LIS_verify(status)
          
! initializing the perturbations to be zero 
          call ESMF_FieldGet(pertField(n),localDE=0,farrayPtr=obs_temp,rc=status)
          call LIS_verify(status)
          obs_temp(:,:) = 0 

          call ESMF_AttributeSet(pertField(n),"Perturbation Type",&
               obs_pert%perttype(1), rc=status)
          call LIS_verify(status)
          
          if(LIS_rc%obs_ngrid(k).gt.0) then 
             call ESMF_AttributeSet(pertField(n),"Standard Deviation",&
                  ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
             call LIS_verify(status)
          endif

          call ESMF_AttributeSet(pertField(n),"Std Normal Max",&
               obs_pert%stdmax(1), rc=status)
          call LIS_verify(status)
          
          call ESMF_AttributeSet(pertField(n),"Ensure Zero Mean",&
               obs_pert%zeromean(1),rc=status)
          call LIS_verify(status)
          
          call ESMF_AttributeSet(pertField(n),"Temporal Correlation Scale",&
               obs_pert%tcorr(1),rc=status)
          
          call ESMF_AttributeSet(pertField(n),"X Correlation Scale",&
               obs_pert%xcorr(1),rc=status)
          
          call ESMF_AttributeSet(pertField(n),"Y Correlation Scale",&
               obs_pert%ycorr(1),rc=status)

          call ESMF_AttributeSet(pertField(n),"Cross Correlation Strength",&
               obs_pert%ccorr(1,:),itemCount=1,rc=status)

          call ESMF_StateAdd(OBS_Pert_State(n),(/pertField(n)/),rc=status)
          call LIS_verify(status)
       endif
          
       deallocate(vname)
       deallocate(varmax)
       deallocate(varmin)
       deallocate(ssdev)
    enddo
!-------------------------------------------------------------
! set up the S1 domain %and interpolation weights. 
!-------------------------------------------------------------

    ! Open first file to get nx ny dimensions
    call system('ls ./' // trim(S1sigmaobsdir) // ' > ./S1_listfiles.txt')

    open(flist, file=trim('./S1_listfiles.txt'), &
         status='old', iostat=status)

    read(flist, '(a)', iostat=status) S1_firstfile
    S1_filename = trim(S1sigmaobsdir) // '/' // trim(S1_firstfile)

    ios = nf90_open(path=S1_filename,&
                    mode=NF90_NOWRITE,ncid=ncid)
    call LIS_verify(ios,&
        'Error reading in S1 data dimensions: Error opening file '// S1_filename)
    ios = nf90_inquire_dimension(ncid,1,yname,NY)
    ios = nf90_inquire_dimension(ncid,2,xname,NX)
    ios = nf90_close(ncid)

    do n=1,LIS_rc%nnest
       S1_sigma_struc(n)%nc = NX
       S1_sigma_struc(n)%nr = NY

       allocate(S1_sigma_struc(n)%sigma(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))
       allocate(S1_sigma_struc(n)%sigmatime(&
            LIS_rc%obs_lnc(k), LIS_rc%obs_lnr(k)))
       S1_sigma_struc(n)%sigma = LIS_rc%udef
       S1_sigma_struc(n)%sigmatime = -1

       call LIS_registerAlarm("S1 backscatter read alarm",&
            86400.0, 86400.0)

       S1_sigma_struc(n)%startMode = .true. 
    enddo

    write(LIS_logunit,*) '[INFO] Created ESMF States to hold S1 observations data'

  end subroutine S1_sigmaVVSM_setup
  
end module S1_sigmaVVSM_Mod
