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
! !MODULE: ANSASCFsnow_Mod
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle ANSA (AFWA NASA snow algorithm) SCF retrievals. 
!   The SCF data is primarily based on AMSR-E, with observation
!   gaps filled in using previous overpasses. For more details, please see 
!   the ANSA reference: 
!
!   Foster et al. ``A blended global snow product using visible, 
!   passive microwave and scatterometer satellite data'', International
!   Journal of Remote Sensing, 2010. 
!   
! !REVISION HISTORY: 
!  1 Jun 09   Sujay Kumar;   Initial Specification
!  1  May 2014 Yuqiong Liu;  adapted for ANSA SCF assimilation using either EnKF or DI
 
module ANSASCFsnow_Mod
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
!EOP
  implicit none
  
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: ANSASCFsnow_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: ANSASCFsnow_struc

  type, public ::  ANSASCFsnow_dec
     logical             :: startMode
     real                :: ssdev
     integer             :: mi
     real                :: cornerlat1,cornerlat2
     real                :: cornerlon1,cornerlon2
     real                :: minlat,minlon
     real                :: maxlat,maxlon
     integer             :: nc,nr
     integer             :: offset1, offset2
     real                :: gridDesc(6)

     integer, allocatable    :: n11(:)
     integer, allocatable    :: n12(:)
     integer, allocatable    :: n21(:)
     integer, allocatable    :: n22(:)
     real,    allocatable    :: w11(:)
     real,    allocatable    :: w12(:)
     real,    allocatable    :: w21(:)
     real,    allocatable    :: w22(:)

     integer, allocatable    :: n113(:)

     real,    allocatable    :: sca(:)
     real                :: assim_lhour
     character(50)       :: scaname
     character(100)      :: scf_fn_conv
     integer             :: obserr
     integer             :: useEnKFwithDI
     character(20)       :: di_opt
     real                :: di_addswe
     real                :: di_maxmelt
     real                :: di_minswe
     real                :: di_ndaymelt
     real                :: di_scaobs1, di_scaobs2, di_scaobs3
     real                :: di_swemod

  end type ANSASCFsnow_dec

  type(ANSASCFsnow_dec), allocatable :: ANSASCFsnow_struc(:)

contains
!BOP
! 
! !ROUTINE: ANSASCFsnow_setup
! \label{ANSASCFsnow_setup}
! 
! !INTERFACE: 
  subroutine ANSASCFsnow_setup(k, OBS_State, OBS_Pert_State)
! !USES: 

    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_historyMod
    use LIS_perturbMod
    use LIS_DAobservationsMod
    use LIS_logMod
   
    implicit none 

! !ARGUMENTS: 
    integer                ::  k 
    type(ESMF_State)       ::  OBS_State(LIS_rc%nnest)
    type(ESMF_State)       ::  OBS_Pert_State(LIS_rc%nnest)
! 
! !DESCRIPTION: 
!   
!   This routine completes the runtime initializations and 
!   creation of data strctures required for ANSA SCF data. 
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
    character(len=LIS_CONST_PATH_LEN) ::  ansasnowobsdir
    character*100          ::  temp
    real,  allocatable         ::  obsstd(:)
    character*1            ::  vid(2)
    character*40, allocatable  ::  vname(:)
    real,         allocatable  ::  varmin(:)
    real,         allocatable  ::  varmax(:)
    type(pert_dec_type)    :: obs_pert
    real, pointer          :: obs_temp(:,:)
    real                   :: gridDesci(LIS_rc%nnest,50)
    real, allocatable          :: ssdev(:)
    real                   :: dx, dy

    allocate(ANSASCFsnow_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"ANSA SCF data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ansasnowobsdir,&
            rc=status)
       call LIS_verify(status,'ANSA snow data directory: not defined')
       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            ansasnowobsdir, rc=status)
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

    write(LIS_logunit,*)'read ANSA SCF data specifications'

!----------------------------------------------------------------------------
!   Create the array containers that will contain the observations and
!   the perturbations. For this synthetic case, it is assumed that the 
!   observations are in the grid space. Since there is only one layer
!   being assimilated, the array size is LIS_rc%ngrid(n). 
!   
!----------------------------------------------------------------------------

    do n=1,LIS_rc%nnest
       
       write(unit=temp,fmt='(i2.2)') 1
       read(unit=temp,fmt='(2a1)') vid

       obsField(n) = ESMF_FieldCreate(grid=LIS_obsVecGrid(n,k),&
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
          ANSASCFsnow_struc(n)%ssdev =obs_pert%ssdev(1) 

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
          call LIS_verify(status)
          
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
! set up the ANSA domain and interpolation weights. 
!-------------------------------------------------------------
    gridDesci = 0 
    call ESMF_ConfigFindLabel(LIS_config,"ANSA SCF lower left lat:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ANSASCFsnow_struc(n)%gridDesc(1),&
            rc=status)
       call LIS_verify(status,'ANSA SCF lower left lat: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ANSA SCF lower left lon:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ANSASCFsnow_struc(n)%gridDesc(2),&
            rc=status)
       call LIS_verify(status,'ANSA SCF lower left lon: not defined')
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config,"ANSA SCF upper right lat:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ANSASCFsnow_struc(n)%gridDesc(3),&
            rc=status)
       call LIS_verify(status,'ANSA SCF upper right lat: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ANSA SCF upper right lon:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ANSASCFsnow_struc(n)%gridDesc(4),&
            rc=status)
       call LIS_verify(status,'ANSA SCF upper right lon: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ANSA SCF resolution (dx):",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ANSASCFsnow_struc(n)%gridDesc(5),&
            rc=status)
       call LIS_verify(status,'ANSA SCF resolution (dx): not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ANSA SCF resolution (dy):",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ANSASCFsnow_struc(n)%gridDesc(6),&
            rc=status)
       call LIS_verify(status,'ANSA SCF resolution (dy): not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ANSA SCF local time for assimilation:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ANSASCFsnow_struc(n)%assim_lhour,&
            rc=status)
       call LIS_verify(status,'ANSA SCF local time for assimilation: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ANSA SCF field name:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ANSASCFsnow_struc(n)%scaname,&
            rc=status)
       call LIS_verify(status,'ANSA SCF field name: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ANSA SCF file name convention:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ANSASCFsnow_struc(n)%scf_fn_conv,&
            rc=status)
       call LIS_verify(status,'ANSA SCF file name convention: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ANSA SCF use triangular-shaped observation error:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ANSASCFsnow_struc(n)%obserr,&
            rc=status)
       call LIS_verify(status,'ANSA SCF use triangular-shaped observation error: not defined')
    enddo 

    call ESMF_ConfigFindLabel(LIS_config,"ANSA SCF using EnKF with DI:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetattribute(LIS_config,ANSASCFsnow_struc(n)%useEnKFwithDI,&
            rc=status)
       call LIS_verify(status,'ANSA SCF using EnKF with DI: not defined')
    enddo   

    call ESMF_ConfigFindLabel(LIS_config,"ANSA SCF direct insertion methodology:",&
         rc=status)
    do n=1,LIS_rc%nnest
       if (ANSASCFsnow_struc(n)%useEnKFwithDI .eq. 1) then
          call ESMF_ConfigGetattribute(LIS_config,ANSASCFsnow_struc(n)%di_opt,&
               rc=status)
          call LIS_verify(status,'ANSA SCF direct insertion methodology: not defined')
       endif
    enddo        

    call ESMF_ConfigFindLabel(LIS_config,"ANSA SCF amount of SWE (mm) to add to model:",&
         rc=status)
    do n=1,LIS_rc%nnest
       if (ANSASCFsnow_struc(n)%useEnKFwithDI .eq. 1) then
          call ESMF_ConfigGetattribute(LIS_config,ANSASCFsnow_struc(n)%di_addswe,&
               rc=status)
          call LIS_verify(status,'ANSA SCF amount of SWE (mm) to add to model: not defined')
       endif
    enddo      

    call ESMF_ConfigFindLabel(LIS_config,"ANSA SCF maximum SWE melt rate (mm/day):",&
         rc=status)
    do n=1,LIS_rc%nnest
       if (ANSASCFsnow_struc(n)%useEnKFwithDI .eq. 1) then
          call ESMF_ConfigGetattribute(LIS_config,ANSASCFsnow_struc(n)%di_maxmelt,&
               rc=status)
          call LIS_verify(status,'ANSA SCF maximum SWE melt rate (mm/day): not defined')
       endif
    enddo    

    call ESMF_ConfigFindLabel(LIS_config,"ANSA SCF threshold of model SWE to be removed at once:",&
         rc=status)
    do n=1,LIS_rc%nnest
       if (ANSASCFsnow_struc(n)%useEnKFwithDI .eq. 1) then
          call ESMF_ConfigGetattribute(LIS_config,ANSASCFsnow_struc(n)%di_minswe,&
               rc=status)
          call LIS_verify(status,'ANSA SCF threshold of model SWE to be removed at once: not defined')
       endif
    enddo  

    call ESMF_ConfigFindLabel(LIS_config,"ANSA SCF length of snowmelt period in days:",&
         rc=status)
    do n=1,LIS_rc%nnest
       if (ANSASCFsnow_struc(n)%useEnKFwithDI .eq. 1) then
          call ESMF_ConfigGetattribute(LIS_config,ANSASCFsnow_struc(n)%di_ndaymelt,&
               rc=status)
          call LIS_verify(status,'ANSA SCF length of snowmelt period in days: not defined')
       endif
    enddo 

    call ESMF_ConfigFindLabel(LIS_config,"ANSA SCF threshold of observed SCF for snow presence:",&
         rc=status)
    do n=1,LIS_rc%nnest
       if (ANSASCFsnow_struc(n)%useEnKFwithDI .eq. 1) then
          call ESMF_ConfigGetattribute(LIS_config,ANSASCFsnow_struc(n)%di_scaobs1,&
               rc=status)
          call LIS_verify(status,'ANSA SCF threshold of observed SCF for snow presence: not defined')
       endif
    enddo 

    call ESMF_ConfigFindLabel(LIS_config,"ANSA SCF threshold of observed SCF for snow non-presence:",&
         rc=status)
    do n=1,LIS_rc%nnest
       if (ANSASCFsnow_struc(n)%useEnKFwithDI .eq. 1) then
          call ESMF_ConfigGetattribute(LIS_config,ANSASCFsnow_struc(n)%di_scaobs2,&
               rc=status)
          call LIS_verify(status,'ANSA SCF threshold of observed SCF for snow non-presence: not defined')
       endif
    enddo 

    call ESMF_ConfigFindLabel(LIS_config,"ANSA SCF threshold of observed SCF for non-full snow cover:",&
         rc=status)
    do n=1,LIS_rc%nnest
       if (ANSASCFsnow_struc(n)%useEnKFwithDI .eq. 1) then
          call ESMF_ConfigGetattribute(LIS_config,ANSASCFsnow_struc(n)%di_scaobs3,&
               rc=status)
          call LIS_verify(status,'ANSA SCF threshold of observed SCF for non-full snow cover: not defined')
       endif
    enddo 

    call ESMF_ConfigFindLabel(LIS_config,"ANSA SCF threshold of model SWE(mm) for snow non-presence:",&
         rc=status)
    do n=1,LIS_rc%nnest
       if (ANSASCFsnow_struc(n)%useEnKFwithDI .eq. 1) then
          call ESMF_ConfigGetattribute(LIS_config,ANSASCFsnow_struc(n)%di_swemod,&
               rc=status)
          call LIS_verify(status,'ANSA SCF threshold of model SWE(mm) for snow non-presence: not defined')
       endif
    enddo 

    do n=1,LIS_rc%nnest

       ANSASCFsnow_struc(n)%minlat = ANSASCFsnow_struc(n)%gridDesc(1)
       ANSASCFsnow_struc(n)%minlon = ANSASCFsnow_struc(n)%gridDesc(2)
       ANSASCFsnow_struc(n)%maxlat = ANSASCFsnow_struc(n)%gridDesc(3)
       ANSASCFsnow_struc(n)%maxlon = ANSASCFsnow_struc(n)%gridDesc(4)
       dx = ANSASCFsnow_struc(n)%gridDesc(5)
       dy = ANSASCFsnow_struc(n)%gridDesc(6)
    enddo
    
    do n=1, LIS_rc%nnest
!sets the local domain corner points with additional buffer 
       ANSASCFsnow_struc(n)%cornerlat1 = max(ANSASCFsnow_struc(n)%minlat, &
            nint((LIS_domain(n)%minlat-ANSASCFsnow_struc(n)%minlat)/dx)*dx+ANSASCFsnow_struc(n)%minlat-2*dx)
       ANSASCFsnow_struc(n)%cornerlon1 = max(ANSASCFsnow_struc(n)%minlon, &
            nint((LIS_domain(n)%minlon-ANSASCFsnow_struc(n)%minlon)/dy)*dy+ANSASCFsnow_struc(n)%minlon-2*dy)
       ANSASCFsnow_struc(n)%cornerlat2 = min(ANSASCFsnow_struc(n)%maxlat, &
            nint((LIS_domain(n)%maxlat-ANSASCFsnow_struc(n)%minlat)/dx)*dx+ANSASCFsnow_struc(n)%minlat+2*dx)
       ANSASCFsnow_struc(n)%cornerlon2 = min(ANSASCFsnow_struc(n)%maxlon, &
            nint((LIS_domain(n)%maxlon-ANSASCFsnow_struc(n)%minlon)/dy)*dy+ANSASCFsnow_struc(n)%minlon+2*dy)


       ANSASCFsnow_struc(n)%offset1 = &
            nint((ANSASCFsnow_struc(n)%cornerlon1-&
            ANSASCFsnow_struc(n)%minlon)/dy)
       ANSASCFsnow_struc(n)%offset2 = &
            nint((ANSASCFsnow_struc(n)%cornerlat1-&
            ANSASCFsnow_struc(n)%minlat)/dx)

       ANSASCFsnow_struc(n)%nr = &
            nint((ANSASCFsnow_struc(n)%cornerlat2-&
            ANSASCFsnow_struc(n)%cornerlat1)/dx)+1
       ANSASCFsnow_struc(n)%nc = &
            nint((ANSASCFsnow_struc(n)%cornerlon2-&
            ANSASCFsnow_struc(n)%cornerlon1)/dy)+1

       gridDesci(n,1) = 0 
       gridDesci(n,2) = ANSASCFsnow_struc(n)%nc
       gridDesci(n,3) = ANSASCFsnow_struc(n)%nr
       gridDesci(n,4) = ANSASCFsnow_struc(n)%cornerlat1
       gridDesci(n,5) = ANSASCFsnow_struc(n)%cornerlon1
       gridDesci(n,6) = 128
       gridDesci(n,7) = ANSASCFsnow_struc(n)%cornerlat2
       gridDesci(n,8) = ANSASCFsnow_struc(n)%cornerlon2
       gridDesci(n,9) = ANSASCFsnow_struc(n)%gridDesc(5)
       gridDesci(n,10) = ANSASCFsnow_struc(n)%gridDesc(6)
       gridDesci(n,20) = 64.0

       ANSASCFsnow_struc(n)%mi = ANSASCFsnow_struc(n)%nc*ANSASCFsnow_struc(n)%nr

       allocate(ANSASCFsnow_struc(n)%n11(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(ANSASCFsnow_struc(n)%n12(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(ANSASCFsnow_struc(n)%n21(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(ANSASCFsnow_struc(n)%n22(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

       allocate(ANSASCFsnow_struc(n)%w11(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(ANSASCFsnow_struc(n)%w12(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(ANSASCFsnow_struc(n)%w21(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(ANSASCFsnow_struc(n)%w22(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

       allocate(ANSASCFsnow_struc(n)%sca(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

       ANSASCFsnow_struc(n)%sca = LIS_rc%udef

       call bilinear_interp_input(n, gridDesci(n,:), &
            ANSASCFsnow_struc(n)%n11,ANSASCFsnow_struc(n)%n12,&
            ANSASCFsnow_struc(n)%n21,ANSASCFsnow_struc(n)%n22,&
            ANSASCFsnow_struc(n)%w11,ANSASCFsnow_struc(n)%w12,&
            ANSASCFsnow_struc(n)%w21,ANSASCFsnow_struc(n)%w22)

       allocate(ANSASCFsnow_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

#if 0 
       call neighbor_interp_input(gridDesci(n,:),&
            LIS_rc%lnc(n), LIS_rc%lnr(n),&
            LIS_domain(n)%lat, LIS_domain(n)%lon,&
            ANSASCFsnow_struc(n)%n113)

#endif
       call LIS_registerAlarm("ANSA SCF read alarm",&
            86400.0, 86400.0)

       ANSASCFsnow_struc(n)%startMode = .true.  
 
    enddo

    write(LIS_logunit,*) 'Created ESMF States to hold ANSA observations data'

  end subroutine ANSASCFsnow_setup
  
end module ANSASCFsnow_Mod

