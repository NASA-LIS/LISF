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
! !MODULE: GLASSalbedo_Mod
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle 
! 
! !REVISION HISTORY: 
!  21 Dec 2017    Sujay Kumar; initial specification
! 
module GLASSalbedo_Mod
! !USES: 
  use ESMF
  use map_utils
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: GLASSalbedo_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: GLASSalbedo_struc
!EOP
  type, public:: GLASSalbedo_dec

     character*100          :: source
     logical                :: startMode
     integer                :: nc
     integer                :: nr
     integer                :: mi
     real*8                 :: time1, time2
     integer                :: fnd
     integer                :: qcflag
     real,     allocatable  :: obs_bs1(:)
     real,     allocatable  :: obs_bs2(:)
     real,     allocatable  :: obs_ws1(:)
     real,     allocatable  :: obs_ws2(:)
     real                   :: gridDesci(50)    
     
     real,    allocatable :: rlat(:)
     real,    allocatable :: rlon(:)
     integer, allocatable :: n11(:)
     integer, allocatable :: n12(:)
     integer, allocatable :: n21(:)
     integer, allocatable :: n22(:)
     real,    allocatable :: w11(:)
     real,    allocatable :: w12(:)
     real,    allocatable :: w21(:)
     real,    allocatable :: w22(:)

  end type GLASSalbedo_dec
  
  type(GLASSalbedo_dec),allocatable :: GLASSalbedo_struc(:)
  
contains

!BOP
! 
! !ROUTINE: GLASSalbedo_setup
! \label{GLASSalbedo_setup}
! 
! !INTERFACE: 
  subroutine GLASSalbedo_setup(k, OBS_State, OBS_Pert_State)
! !USES: 
    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_historyMod
    use LIS_dataAssimMod
    use LIS_perturbMod
    use LIS_DAobservationsMod
    use LIS_logmod

    implicit none 

! !ARGUMENTS: 
    integer                ::  k
    type(ESMF_State)       ::  OBS_State(LIS_rc%nnest)
    type(ESMF_State)       ::  OBS_Pert_State(LIS_rc%nnest)
! 
! !DESCRIPTION: 
!   
!   This routine completes the runtime initializations and 
!   creation of data strctures required for handling NASASMAP 
!   soil moisture data. 
!  
!   The arguments are: 
!   \begin{description}
!    \item[OBS\_State]   observation state 
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP
    integer                ::  n,m,i,t,kk,jj
    integer                ::  ftn
    integer                ::  status
    type(ESMF_Field)       ::  obsField
    type(ESMF_Field)       ::  pertField
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_ArraySpec)   ::  pertArrSpec
    character(len=LIS_CONST_PATH_LEN) ::  albedoobsdir
    character*100          ::  temp
    real,  allocatable         ::  ssdev(:)
    real,  allocatable         ::  obsstd(:)
    character*1            ::  vid(2)
    type(pert_dec_type)    ::  obs_pert
    real, pointer          ::  obs_temp(:,:)
    character*40, allocatable  ::  vname(:)
    real        , allocatable  ::  varmin(:)
    real        , allocatable  ::  varmax(:)
    integer                :: c,r
    real                   ::  cornerlat1, cornerlon1
    real                   ::  cornerlat2, cornerlon2

    allocate(GLASSalbedo_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"GLASS Albedo data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,albedoobsdir,&
            rc=status)
       call LIS_verify(status, 'GLASS Albedo data directory: is missing')

       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            albedoobsdir, rc=status)
       call LIS_verify(status)
    enddo

! source = "AVHRR" or "MODIS"
    call ESMF_ConfigFindLabel(LIS_config,"GLASS Albedo data source:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,GLASSalbedo_struc(n)%source,&
            rc=status)
       call LIS_verify(status, 'GLASS Albedo data source: is missing')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"GLASS Albedo apply QC flags:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,GLASSalbedo_struc(n)%qcflag,&
            rc=status)
       call LIS_verify(status, 'GLASS Albedo apply QC flags: is missing')
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

    write(LIS_logunit,*)&
         '[INFO] read GLASS Albedo data specifications'       

    do n=1,LIS_rc%nnest

       allocate(GLASSalbedo_struc(n)%obs_bs1(&
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
       allocate(GLASSalbedo_struc(n)%obs_bs2(&
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
       allocate(GLASSalbedo_struc(n)%obs_ws1(&
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
       allocate(GLASSalbedo_struc(n)%obs_ws2(&
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
       
       do m=1,LIS_rc%nobtypes(k)
          write(unit=temp,fmt='(i2.2)') m
          read(unit=temp,fmt='(2a1)') vid
          
          obsField = ESMF_FieldCreate(arrayspec=realarrspec,&
               grid=LIS_obsVecGrid(n,k),&
               name="Observation"//vid(1)//vid(2), rc=status)
          call LIS_verify(status)

          call ESMF_StateAdd(OBS_State(n),(/obsField/),rc=status)
          call LIS_verify(status)
       enddo

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

       allocate(ssdev(LIS_rc%obs_ngrid(k)))
      
       if(trim(LIS_rc%perturb_obs(k)).ne."none") then 
          allocate(obs_pert%vname(LIS_rc%nobtypes(k)))
          allocate(obs_pert%perttype(LIS_rc%nobtypes(k)))
          allocate(obs_pert%ssdev(LIS_rc%nobtypes(k)))
          allocate(obs_pert%stdmax(LIS_rc%nobtypes(k)))
          allocate(obs_pert%zeromean(LIS_rc%nobtypes(k)))
          allocate(obs_pert%tcorr(LIS_rc%nobtypes(k)))
          allocate(obs_pert%xcorr(LIS_rc%nobtypes(k)))
          allocate(obs_pert%ycorr(LIS_rc%nobtypes(k)))
          allocate(obs_pert%ccorr(LIS_rc%nobtypes(k),&
               LIS_rc%nobtypes(k)))
          
          call LIS_readPertAttributes(LIS_rc%nobtypes(k),&
               LIS_rc%obspertAttribfile(k),&
               obs_pert)
          
          ! Set obs err to be uniform (will be rescaled later for each grid point). 
          do m=1,LIS_rc%nobtypes(k)
             ssdev = obs_pert%ssdev(m)
          
             write(unit=temp,fmt='(i2.2)') m
             read(unit=temp,fmt='(2a1)') vid
             pertField = ESMF_FieldCreate(arrayspec=pertArrSpec,&
                  grid=LIS_obsEnsOnGrid(n,k),&
                  name="Observation"//vid(1)//vid(2),&
                  rc=status)
             call LIS_verify(status)
          
             ! initializing the perturbations to be zero 
             call ESMF_FieldGet(pertField,localDE=0,&
                  farrayPtr=obs_temp,rc=status)
             call LIS_verify(status)
             obs_temp(:,:) = 0 
             
             call ESMF_AttributeSet(pertField,"Perturbation Type",&
                  obs_pert%perttype(m), rc=status)
             call LIS_verify(status)
             
             if(LIS_rc%obs_ngrid(k).gt.0) then 
                call ESMF_AttributeSet(pertField,"Standard Deviation",&
                     ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
                call LIS_verify(status)
             endif
             
             call ESMF_AttributeSet(pertField,"Std Normal Max",&
                  obs_pert%stdmax(m), rc=status)
             call LIS_verify(status)
             
             call ESMF_AttributeSet(pertField,"Ensure Zero Mean",&
                  obs_pert%zeromean(m),rc=status)
             call LIS_verify(status)
             
             call ESMF_AttributeSet(pertField,"Temporal Correlation Scale",&
                  obs_pert%tcorr(m),rc=status)
             call LIS_verify(status)
             
             call ESMF_AttributeSet(pertField,"X Correlation Scale",&
                  obs_pert%xcorr(m),rc=status)
             
             call ESMF_AttributeSet(pertField,"Y Correlation Scale",&
                  obs_pert%ycorr(m),rc=status)
             
             call ESMF_AttributeSet(pertField,"Cross Correlation Strength",&
                  obs_pert%ccorr(m,:),itemCount=1,rc=status)
             
             call ESMF_StateAdd(OBS_Pert_State(n),(/pertField/),rc=status)
             call LIS_verify(status)         
          enddo
       endif
       
       deallocate(vname)
       deallocate(varmax)
       deallocate(varmin)
       deallocate(ssdev)   
       
    enddo
    
    do n=1,LIS_rc%nnest
          
       cornerlat1 = max(-59.975, nint((LIS_rc%obs_gridDesc(k,4)+59.975)/0.05)*0.05-59.975-2*0.05)
       cornerlon1 = max(-179.975, nint((LIS_rc%obs_gridDesc(k,5)+179.975)/0.05)*0.05-179.975-2*0.05)
       cornerlat2 = min(89.975, nint((LIS_rc%obs_gridDesc(k,7)+59.975)/0.05)*0.05-59.975+2*0.05)
       cornerlon2 = min(179.975, nint((LIS_rc%obs_gridDesc(k,8)+179.975)/0.05)*0.05-179.975+2*0.05)


       GLASSalbedo_struc(n)%nr = nint((cornerlat2-cornerlat1)/0.05)+1
       GLASSalbedo_struc(n)%nc = nint((cornerlon2-cornerlon1)/0.05)+1

       GLASSalbedo_struc(n)%gridDesci(1) = 0 
       GLASSalbedo_struc(n)%gridDesci(2) = GLASSalbedo_struc(n)%nc
       GLASSalbedo_struc(n)%gridDesci(3) = GLASSalbedo_struc(n)%nr 
       GLASSalbedo_struc(n)%gridDesci(4) = cornerlat1
       GLASSalbedo_struc(n)%gridDesci(5) = cornerlon1
       GLASSalbedo_struc(n)%gridDesci(6) = 128
       GLASSalbedo_struc(n)%gridDesci(7) = cornerlat2
       GLASSalbedo_struc(n)%gridDesci(8) = cornerlon2
       GLASSalbedo_struc(n)%gridDesci(9) = 0.05
       GLASSalbedo_struc(n)%gridDesci(10) = 0.05
       GLASSalbedo_struc(n)%gridDesci(20) = 64

       GLASSalbedo_struc(n)%mi = GLASSalbedo_struc(n)%nc*GLASSalbedo_struc(n)%nr

!-----------------------------------------------------------------------------
!   Use interpolation if LIS is running finer than 5km. 
!-----------------------------------------------------------------------------
       if(LIS_rc%obs_gridDesc(k,10).le.0.05) then 

          allocate(GLASSalbedo_struc(n)%rlat(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(GLASSalbedo_struc(n)%rlon(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))

          allocate(GLASSalbedo_struc(n)%n11(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(GLASSalbedo_struc(n)%n12(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(GLASSalbedo_struc(n)%n21(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(GLASSalbedo_struc(n)%n22(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(GLASSalbedo_struc(n)%w11(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(GLASSalbedo_struc(n)%w12(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(GLASSalbedo_struc(n)%w21(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(GLASSalbedo_struc(n)%w22(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          
          call bilinear_interp_input_withgrid(&
               GLASSalbedo_struc(n)%gridDesci(:), &
               LIS_rc%obs_gridDesc(k,:),&
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),&
               GLASSalbedo_struc(n)%rlat, GLASSalbedo_struc(n)%rlon,&
               GLASSalbedo_struc(n)%n11, GLASSalbedo_struc(n)%n12, &
               GLASSalbedo_struc(n)%n21, GLASSalbedo_struc(n)%n22, &
               GLASSalbedo_struc(n)%w11, GLASSalbedo_struc(n)%w12, &
               GLASSalbedo_struc(n)%w21, GLASSalbedo_struc(n)%w22)
       else
          
          allocate(GLASSalbedo_struc(n)%n11(&
               GLASSalbedo_struc(n)%nc*GLASSalbedo_struc(n)%nr))

          call upscaleByAveraging_input(GLASSalbedo_struc(n)%gridDesci(:),&
               LIS_rc%obs_gridDesc(k,:),&
               GLASSalbedo_struc(n)%nc*GLASSalbedo_struc(n)%nr, &
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), GLASSalbedo_struc(n)%n11)
       endif
       
       call LIS_registerAlarm("GLASS Albedo read alarm",&
            86400.0, 86400.0)

       GLASSalbedo_struc(n)%startMode = .true. 

    enddo
  end subroutine GLASSalbedo_setup
end module GLASSalbedo_Mod
