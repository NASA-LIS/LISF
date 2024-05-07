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
! !MODULE: AmerifluxobsMod
! 
! !DESCRIPTION: 
!
!   
! !REVISION HISTORY: 
!  09 Jul 09    Sujay Kumar;   Initial Specification
! 
module AmerifluxobsMod
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
!EOP
  implicit none
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: Amerifluxobs_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: Amerifluxobs_struc

  type, public ::  Amerifluxobs_data_dec
     character(len=LIS_CONST_PATH_LEN) :: odir
     integer                 :: n_stns
     character*100, allocatable  :: site_name(:)
     character*100, allocatable  :: stn_name(:)
     real,          allocatable  :: stnlat(:)
     real,          allocatable  :: stnlon(:)
     real,          allocatable  :: sfsm_wt(:,:)
     real,          allocatable  :: rzsm_wt(:,:)
     real,          allocatable  :: sfst_wt(:,:)
     real,          allocatable  :: rzst_wt(:,:)
     logical                 :: startflag
     integer                 :: yr
     integer                 :: nsmlayers
     integer                 :: nstlayers

     type(ESMF_Time)         :: starttime
     type(ESMF_TimeInterval) :: timestep

     real,          allocatable  :: Qle(:,:)        ! Latent heat flux
     real,          allocatable  :: Qh(:,:)         ! Sensible heat flux
     real,          allocatable  :: Qg(:,:)         ! Ground Heat Flux
     real,          allocatable  :: sfsm(:,:)
     real,          allocatable  :: sfst(:,:)
  end type Amerifluxobs_data_dec

  type(Amerifluxobs_data_dec), allocatable :: Amerifluxobs_struc(:)

contains
!BOP
! 
! !ROUTINE: Amerifluxobs_setup
! \label{Amerifluxobs_setup}
! 
! !INTERFACE: 
  subroutine Amerifluxobs_setup(Obs_State)
! !USES: 
    use LIS_coreMod, only : LIS_rc, LIS_config, LIS_vecGrid
    use LIS_logMod, only : LIS_logunit, LIS_verify, & 
         LIS_getNextUnitNumber, LIS_releaseUnitNumber

    implicit none 

! !ARGUMENTS: 
    type(ESMF_State)       ::  Obs_State(LIS_rc%nnest)
! 
! !DESCRIPTION: 
!   
!   The arguments are: 
!   \begin{description}
!    \item[Obj\_Space]   observation/Objective space object 
!   \end{description}
!EOP
    integer                   ::  n 
    integer                   ::  status
    integer                   ::  i 
    type(ESMF_ArraySpec)      ::  realarrspec
    type(ESMF_Field)          ::  obsField
    character(len=LIS_CONST_PATH_LEN) ::  obsdir
    character*100             ::  vname
    character*200             ::  currentLine
    real                      ::  swcd(2), tsd(2)
    integer                   ::  k,iloc
    integer                   ::  arrayLen
    character(len=LIS_CONST_PATH_LEN) :: obsAttribFile(LIS_rc%nnest)
    integer                   ::  ftn
    real                      ::  gridDesci(LIS_rc%nnest,50)
    character(len=LIS_CONST_PATH_LEN) :: stnlist_file

    allocate(Amerifluxobs_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)
    
    call ESMF_ConfigFindLabel(LIS_config,"Ameriflux data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,obsdir,&
            rc=status)
       call LIS_verify(status, 'Ameriflux data directory: not defined')

       call ESMF_AttributeSet(Obs_State(n),"Data Directory",&
            obsdir, rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(Obs_State(n),"Data Update Status",&
            .false., rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config, "Ameriflux station list file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,stnlist_file,&
            rc=status)
       call LIS_verify(status, 'Ameriflux station list file: not defined')
    enddo

    do n=1,LIS_rc%nnest
       ftn = LIS_getNextUnitNumber()
       write(LIS_logunit,*) 'Reading Ameriflux station list file ',&
            trim(stnlist_file)
       open(ftn, file=trim(stnlist_file), form='formatted')
       read(ftn,*)
       read(ftn,*) Amerifluxobs_struc(n)%n_stns
       read(ftn,*) 
       
       allocate(Amerifluxobs_struc(n)%site_name(Amerifluxobs_struc(n)%n_stns))
       allocate(Amerifluxobs_struc(n)%stn_name(Amerifluxobs_struc(n)%n_stns))
       allocate(Amerifluxobs_struc(n)%stnlat(Amerifluxobs_struc(n)%n_stns))
       allocate(Amerifluxobs_struc(n)%stnlon(Amerifluxobs_struc(n)%n_stns))
       allocate(Amerifluxobs_struc(n)%sfsm_wt(Amerifluxobs_struc(n)%n_stns,2))
       allocate(Amerifluxobs_struc(n)%sfst_wt(Amerifluxobs_struc(n)%n_stns,2))

!---------------------------------------------------------------------------
! Leap year test, following rules by Microsoft. The depth of this test may not
! be entirely necessary due to how recent the data is, but it is correct.
!---------------------------------------------------------------------------   
       if (Mod(LIS_rc%yr, 4) == 0) Then            !Step1
          if (Mod(LIS_rc%yr, 100) == 0) Then      !Step2
             if (Mod(LIS_rc%yr, 400) == 0) Then   !Step3
                arrayLen = 17568                  !Step4, leap LIS_rc%yr
             else                           
                arrayLen = 17520                  !Step5, not leap LIS_rc%yr
             end if
          else
             arrayLen = 17568                     !Step4, leap LIS_rc%yr
          end if
       else
          arrayLen = 17520                        !Step5, not leap LIS_rc%yr
       End if
       !svk : to be safe:
       arrayLen = 17570

       allocate(Amerifluxobs_struc(n)%Qle(Amerifluxobs_struc(n)%n_stns, arrayLen))
       allocate(Amerifluxobs_struc(n)%Qh(Amerifluxobs_struc(n)%n_stns, arrayLen))
       allocate(Amerifluxobs_struc(n)%Qg(Amerifluxobs_struc(n)%n_stns, arrayLen))
       allocate(Amerifluxobs_struc(n)%sfsm(Amerifluxobs_struc(n)%n_stns, arrayLen))
       allocate(Amerifluxobs_struc(n)%sfst(Amerifluxobs_struc(n)%n_stns, arrayLen))
       Amerifluxobs_struc(n)%Qle = LIS_rc%udef
       Amerifluxobs_struc(n)%Qh = LIS_rc%udef
       Amerifluxobs_struc(n)%Qg = LIS_rc%udef
       Amerifluxobs_struc(n)%sfsm = LIS_rc%udef
       Amerifluxobs_struc(n)%sfst = LIS_rc%udef

!------------------------------------------------------------------------------
! For each station, this reads the site name, station name, and station
! position in lat, lon coordinates.
!------------------------------------------------------------------------------
       do k=1,AmerifluxObs_struc(n)%n_stns
          read(ftn,'(a)') currentLine
          iloc = Index(currentLine, ";")
          READ(currentLine(1: iloc - 1), *) AmerifluxObs_struc(n)%site_name(k)
          currentLine = currentLine(iloc + 1: Len(currentLine))
          
          iloc = Index(currentLine, ";")
          READ(currentLine(1: iloc - 1), *) AmerifluxObs_struc(n)%stn_name(k)
          currentLine = currentLine(iloc + 1: Len(currentLine))
          
          iloc = Index(currentLine, ";")
          READ(currentLine(1: iloc - 1), *) AmerifluxObs_struc(n)%stnlat(k)
          currentLine = currentLine(iloc + 1: Len(currentLine))
          
          iloc = Index(currentLine, ";")
          READ(currentLine(1: iloc - 1), *) AmerifluxObs_struc(n)%stnlon(k)
          currentLine = currentLine(iloc + 1: Len(currentLine))
          
          iloc = Index(currentLine, ";")
          READ(currentLine(1: iloc - 1), *) SWCD(1)
          currentLine = currentLine(iloc + 1: Len(currentLine))
          
          iloc = Index(currentLine, ";")
          READ(currentLine(1: iloc - 1), *) SWCD(2)
          currentLine = currentLine(iloc + 1: Len(currentLine))
          
          iloc = Index(currentLine, ";") 
          READ(currentLine(1: iloc - 1), *) tsd(1)
          currentLine = currentLine(iloc + 1: Len(currentLine))
          
          READ(currentLine, *) tsd(2)
          
          if(swcd(1).gt.0.and.swcd(2).gt.0) then 
             AmerifluxObs_struc(n)%nsmlayers = 2
          elseif(swcd(1).gt.0.and.swcd(2).le.0) then 
             AmerifluxObs_struc(n)%nsmlayers = 1
          elseif(swcd(1).le.0.and.swcd(2).gt.0) then 
             AmerifluxObs_struc(n)%nsmlayers = 1
          else
             AmerifluxObs_struc(n)%nsmlayers = 0
          endif
          
          if(tsd(1).gt.0.and.tsd(2).gt.0) then 
             AmerifluxObs_struc(n)%nstlayers = 2
          elseif(swcd(1).gt.0.and.swcd(2).le.0) then 
             AmerifluxObs_struc(n)%nstlayers = 1
          elseif(swcd(1).le.0.and.swcd(2).gt.0) then 
             AmerifluxObs_struc(n)%nstlayers = 1
          else
             AmerifluxObs_struc(n)%nstlayers = 0
          endif

       !print*, AmerifluxObs_struc(n)%site_name(k), ' ', AmerifluxObs_struc(n)%stn_name(k), ' ', AmerifluxObs_struc(n)%stnlat(k), ' ', &
!& AmerifluxObs_struc(n)%stnlon(k), ' ', AmerifluxObs_struc(n)%SWC1D(k), ' ', AmerifluxObs_struc(n)%SWC2D(k)
       
!       if(AmerifluxObs_struc(n)%nsmlayers.ne.0) then 
!          call compute_vinterp_weights(LIS_rc%nsmlayers, &
!               AmerifluxObs_struc(n)%nsmlayers,&
!               LIS_rc%lis_sf_d, &
!               LIS_rc%lis_rz_d, &
!               LIS_rc%smthick(:), &
!               swcd(1:AmerifluxObs_struc(n)%nsmlayers)/100.0, &
!               AmerifluxObs_struc(n)%sfsm_wt(k,:), &
!               AmerifluxObs_struc(n)%rzsm_wt(k,:))
!       else
!             AmerifluxObs_struc(n)%sfsm_wt(k,:) = 0.0
!             AmerifluxObs_struc(n)%rzsm_wt(k,:) = 0.0
!       endif

!       if(AmerifluxObs_struc(n)%nstlayers.ne.0) then 
!          call compute_vinterp_weights(LIS_rc%nstlayers, &
!               AmerifluxObs_struc(n)%nstlayers,&
!               LIS_rc%lis_sf_d, &
!               LIS_rc%lis_rz_d, &
!               LIS_rc%stthick(:), &
!               tsd(1:AmerifluxObs_struc(n)%nstlayers)/100.0, &
!               AmerifluxObs_struc(n)%sfst_wt(k,:), &
!               AmerifluxObs_struc(n)%rzst_wt(k,:))
!       else
!             AmerifluxObs_struc(n)%sfst_wt(k,:) = 0.0
!             AmerifluxObs_struc(n)%rzst_wt(k,:) = 0.0
!       endif
          write(LIS_logunit,*) k, AmerifluxObs_struc(n)%site_name(k), &
               AmerifluxObs_struc(n)%stn_name(k), &
               AmerifluxObs_struc(n)%stnlat(k), AmerifluxObs_struc(n)%stnlon(k)
       end do
       call LIS_releaseUnitNumber(ftn)
    enddo

    
    call ESMF_ConfigFindLabel(LIS_config,"Ameriflux observations attributes file:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,obsAttribFile(n),rc=status)
       call LIS_verify(status, 'Ameriflux observations attributes file: not defined')
   
       ftn=LIS_getNextUnitNumber()
       open(ftn,file=trim(obsAttribFile(n)),status='old')
       read(ftn,*)
       
       do k=1,3
          read(ftn,fmt='(a100)') vname
          obsField = ESMF_FieldCreate(arrayspec=realarrspec, grid=LIS_vecGrid(n), &
               name=trim(vname), rc=status)
          call LIS_verify(status)
          
          call ESMF_StateAdd(Obs_State(n),(/obsField/),rc=status)
          call LIS_verify(status)
       enddo

       call LIS_releaseUnitNumber(ftn)
    enddo


    write(LIS_logunit,*) 'Created the States to hold the Ameriflux observations'
    
  end subroutine Amerifluxobs_setup
  
end module AmerifluxobsMod
