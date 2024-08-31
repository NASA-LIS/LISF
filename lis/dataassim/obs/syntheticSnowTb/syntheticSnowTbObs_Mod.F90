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
! !MODULE: syntheticSnowTbObs_Mod
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle synthetic passive microwave brightness temperature (TB) observations. 
!   
! !REVISION HISTORY: 
!  27Feb05    Sujay Kumar;   Initial Specification
!  29Sep17    Yonghwan Kwon;  modified for TB observations
! 
! 
module syntheticSnowTbObs_Mod
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
!EOP
  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: syntheticSnowTbObs_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
contains
!BOP
! 
! !ROUTINE: syntheticSnowTbObs_setup
! \label{syntheticSnowTbObs_setup}
! 
! !INTERFACE: 
  subroutine syntheticSnowTbObs_setup(k, OBS_State, OBS_Pert_State)
! !USES: 
    use LIS_coreMod
    use LIS_logMod
    use LIS_perturbMod
    use LIS_DAobservationsMod
! !ARGUMENTS: 
    integer                ::  k 
    type(ESMF_State)       ::  OBS_State(LIS_rc%nnest)
    type(ESMF_State)       ::  OBS_Pert_State(LIS_rc%nnest)
! 
! !DESCRIPTION:
!   
!   This routine completes the runtime initializations and 
!   creation of data strctures required for TB assimilation
!  
!   The arguments are: 
!   \begin{description}
!    \item[OBS\_State]   observation state 
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP
    integer                ::  n, m 
    integer                ::  ftn
    integer                ::  i
    integer                ::  status
    type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  pertArrSpec
    character(len=LIS_CONST_PATH_LEN) ::  synTB18Vobsdir, synTB18Hobsdir, synTB36Vobsdir, synTB36Hobsdir
    character*100          ::  temp
    real, allocatable          :: ssdev(:)
    character*1            ::  vid(2)
    character*40, allocatable  ::  vname(:)
    real,         allocatable  ::  varmin(:)
    real,         allocatable  ::  varmax(:)
    type(pert_dec_type)    ::  obs_pert
    real, pointer          ::  obs_temp(:,:)

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    !-------------------------------TB_18V-------------------------------------
    call ESMF_ConfigFindLabel(LIS_config,"Synthetic TB_18V data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,synTB18Vobsdir,&
            rc=status)
       call LIS_verify(status)
       call ESMF_AttributeSet(OBS_State(n),"Data Directory TB_18V",&
            synTB18Vobsdir, rc=status)
       call LIS_verify(status)
    enddo
    !-------------------------------TB_18H-------------------------------------
    call ESMF_ConfigFindLabel(LIS_config,"Synthetic TB_18H data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,synTB18Hobsdir,&
            rc=status)
       call LIS_verify(status)
       call ESMF_AttributeSet(OBS_State(n),"Data Directory TB_18H",&
            synTB18Hobsdir, rc=status)
       call LIS_verify(status)
    enddo
    !-------------------------------TB_36V-------------------------------------
    call ESMF_ConfigFindLabel(LIS_config,"Synthetic TB_36V data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,synTB36Vobsdir,&
            rc=status)
       call LIS_verify(status)
       call ESMF_AttributeSet(OBS_State(n),"Data Directory TB_36V",&
            synTB36Vobsdir, rc=status)
       call LIS_verify(status)
    enddo
    !-------------------------------TB_36H-------------------------------------
    call ESMF_ConfigFindLabel(LIS_config,"Synthetic TB_36H data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,synTB36Hobsdir,&
            rc=status)
       call LIS_verify(status)
       call ESMF_AttributeSet(OBS_State(n),"Data Directory TB_36H",&
            synTB36Hobsdir, rc=status)
       call LIS_verify(status)
    enddo
    !------------------------------------------------------------------------
   
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

    write(LIS_logunit,*)'read Synthetic TB data specifications'

!----------------------------------------------------------------------------
!   Create the array containers that will contain the observations and
!   the perturbations. For this synthetic case, it is assumed that the 
!   observations are in the grid space. There are four types of observations
!   (i.e., 18V, 18H, 36V, and 36H) being assimilated.
!   
!----------------------------------------------------------------------------
    do n=1,LIS_rc%nnest
       
       write(unit=temp,fmt='(i2.2)') 1
       read(unit=temp,fmt='(2a1)') vid

       obsField(n) = ESMF_FieldCreate(grid=LIS_obsvecGrid(n,k),&
            arrayspec=realarrspec,&
            name="Observation"//vid(1)//vid(2), rc=status)
       call LIS_verify(status)
       
!Perturbations State
       write(LIS_logunit,*) 'Opening attributes for observations ',&
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
          write(LIS_logunit,*) vname(i),varmin(i),varmax(i)
       enddo
       call LIS_releaseUnitNumber(ftn)  
       
       call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
       call LIS_verify(status)

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
          allocate(obs_pert%ccorr(LIS_rc%nobtypes(k),LIS_rc%nobtypes(k)))

          call LIS_readPertAttributes(LIS_rc%nobtypes(k),&
               LIS_rc%obspertAttribfile(k),&
               obs_pert)

          do m=1,LIS_rc%nobtypes(k)
             write(unit=temp,fmt='(i2.2)') m
             read(unit=temp,fmt='(2a1)') vid
 
             ssdev = obs_pert%ssdev(m)

             pertField(n) = ESMF_FieldCreate(arrayspec=pertArrSpec,&
                  grid=LIS_obsEnsOnGrid(n,k),name="Observation"//vid(1)//vid(2),&
                  rc=status)
             call LIS_verify(status)
          
! initializing the perturbations to be zero 
             call ESMF_FieldGet(pertField(n),localDE=0,farrayPtr=obs_temp,rc=status)
             call LIS_verify(status)
             obs_temp(:,:) = 0 

             call ESMF_AttributeSet(pertField(n),"Perturbation Type",&
                  obs_pert%perttype(m), rc=status)
             call LIS_verify(status)
          
             if(LIS_rc%obs_ngrid(k).gt.0) then 
                call ESMF_AttributeSet(pertField(n),"Standard Deviation",&
                     ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
                call LIS_verify(status)
             endif

             call ESMF_AttributeSet(pertField(n),"Std Normal Max",&
                  obs_pert%stdmax(m), rc=status)
             call LIS_verify(status)
          
             call ESMF_AttributeSet(pertField(n),"Ensure Zero Mean",&
                  obs_pert%zeromean(m),rc=status)
             call LIS_verify(status)
          
             call ESMF_AttributeSet(pertField(n),"Temporal Correlation Scale",&
                  obs_pert%tcorr(m),rc=status)
             call LIS_verify(status)
          
             call ESMF_AttributeSet(pertField(n),"X Correlation Scale",&
                  obs_pert%xcorr(m),rc=status)
          
             call ESMF_AttributeSet(pertField(n),"Y Correlation Scale",&
                  obs_pert%ycorr(m),rc=status)

             call ESMF_AttributeSet(pertField(n),"Cross Correlation Strength",&
                  obs_pert%ccorr(m,:),itemCount=LIS_rc%nobtypes(k),rc=status)

             call ESMF_StateAdd(OBS_Pert_State(n),(/pertField(n)/),rc=status)
             call LIS_verify(status)
          enddo
       endif
          
       deallocate(vname)
       deallocate(varmax)
       deallocate(varmin)
       deallocate(ssdev)
    enddo

    write(LIS_logunit,*) 'Created the States to hold the observations data'
    
  end subroutine syntheticSnowTbObs_setup
  
end module syntheticSnowTbObs_Mod
