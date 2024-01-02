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
! !MODULE: ARMdata_module
! 
! !DESCRIPTION: 
!  
!   
! !REVISION HISTORY: 
!  29 Mar 11    Sujay Kumar;   Initial Specification
! 
module ARMdata_module
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
!EOP
  implicit none
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: ARMdata_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: ARMdata_struc

  type, public ::  ARMdata_data_dec

     character(len=LIS_CONST_PATH_LEN) :: odir
     integer                 :: n_stns 
     character*30            :: site_id
     real                    :: udef 
     real,           allocatable :: stnlat(:)
     real,           allocatable :: stnlon(:)
     character*100,  allocatable :: stn_name(:)
     real,           allocatable :: stnbd(:)

     logical                 :: startFlag

     integer                 :: da
     integer                 :: baebbr_select
     integer                 :: ecor_select
     integer                 :: ebbr_select

     real, allocatable           :: baebbr_qle_p(:,:)
     real, allocatable           :: baebbr_qh_p(:,:)
     real, allocatable           :: baebbr_qg_p(:,:)
     real, allocatable           :: baebbr_qle_c(:,:)
     real, allocatable           :: baebbr_qh_c(:,:)
     real, allocatable           :: baebbr_qg_c(:,:)
     integer, allocatable        :: baebbr_tindex(:,:)

     real, allocatable           :: ebbr_sfsm_p(:,:)
     real, allocatable           :: ebbr_sfst_p(:,:)
     real, allocatable           :: ebbr_sfsm_c(:,:)
     real, allocatable           :: ebbr_sfst_c(:,:)
     integer, allocatable        :: ebbr_tindex(:,:)

     real, allocatable           :: ecor_qle_p(:,:)
     real, allocatable           :: ecor_qh_p(:,:)
     real, allocatable           :: ecor_qle_c(:,:)
     real, allocatable           :: ecor_qh_c(:,:)
     integer, allocatable        :: ecor_tindex(:,:)
     
     type(ESMF_Time)         :: starttime
     type(ESMF_TimeInterval) :: baebbr_ts
     type(ESMF_TimeInterval) :: ecor_ts
     type(ESMF_TimeInterval) :: ebbr_ts

  end type ARMdata_data_dec

  type(ARMdata_data_dec), allocatable :: ARMdata_struc(:)

contains
!BOP
! 
! !ROUTINE: ARMdata_setup
! \label{ARMdata_setup}
! 
! !INTERFACE: 
  subroutine ARMdata_setup(Obj_space)
! !USES: 
    use LIS_coreMod, only : LIS_rc, LIS_config, LIS_vecGrid
    use LIS_logMod, only : LIS_logunit, LIS_verify, & 
         LIS_getNextUnitNumber, LIS_releaseUnitNumber

    implicit none 

! !ARGUMENTS: 
    type(ESMF_State)       ::  Obj_space(LIS_rc%nnest)
! 
! !DESCRIPTION: 
!   
!   This routine completes the runtime initializations and 
!   creation of data strctures required for reading the PBMR
!   soil moisture data for the Walnut Gulch domain. 
!  
!   The arguments are: 
!   \begin{description}
!    \item[Obj\_Space]   observation/Objective space object 
!   \end{description}
!EOP
    integer                   ::  n 
    integer                   ::  status
    integer                   ::  nobstypes
    integer                   ::  i 
    type(ESMF_ArraySpec)      ::  realarrspec
    type(ESMF_Field),allocatable  ::  obsField(:)
    character(len=LIS_CONST_PATH_LEN) ::  smobsdir
    character*100, allocatable    ::  vname(:)
    character(len=LIS_CONST_PATH_LEN) ::  obsAttribFile(LIS_rc%nnest)
    integer                   ::  ftn
    integer                   ::  k, iloc
    character*100             ::  currentLine
    character(len=LIS_CONST_PATH_LEN) ::  stnlist_file
    character(len=LIS_CONST_PATH_LEN) ::  objspaceAttribFile(LIS_rc%nnest)


    allocate(ARMdata_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)
    
    call ESMF_ConfigFindLabel(LIS_config,"ARM data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,ARMdata_struc(n)%odir,&
            rc=status)
       call LIS_verify(status, "ARM data directory: not defined")
       
       call ESMF_AttributeSet(Obj_Space(n),"Data Update Status",&
            .false., rc=status)
       call LIS_verify(status)
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config,"ARM site identifier name:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,ARMdata_struc(n)%site_id,&
            rc=status)
       call LIS_verify(status, "ARM site identifier name: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config, "ARM station list file:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config, stnlist_file,&
            rc=status)
       call LIS_verify(status, "ARM station list file: not defined")
       ftn = LIS_getNextUnitNumber()
       open(ftn,file=trim(stnlist_file), form='formatted')
       read(ftn,*)
       read(ftn,*) ARMdata_struc(n)%n_stns
       read(ftn,*) 
       
       allocate(ARMdata_struc(n)%stn_name(ARMdata_struc(n)%n_stns))
       allocate(ARMdata_struc(n)%stnlat(ARMdata_struc(n)%n_stns))
       allocate(ARMdata_struc(n)%stnlon(ARMdata_struc(n)%n_stns))
       allocate(ARMdata_struc(n)%stnbd(ARMdata_struc(n)%n_stns))

       do k=1,ARMdata_struc(n)%n_stns
          read(ftn,'(a)') currentLine
          iloc = Index(currentLine, ";")
          read(currentLine(1:iloc -1), *) ARMdata_struc(n)%stn_name(k)
          currentLine = currentLine(iloc+1:Len(currentLine))
          
          iloc = Index(currentLine, ";")
          read(currentLine(1:iloc -1), *) ARMdata_struc(n)%stnlat(k)
          currentLine = currentLine(iloc+1:Len(currentLine))

          iloc = Index(currentLine, ";")
          read(currentLine(1:iloc -1), *) ARMdata_struc(n)%stnlon(k)
          currentLine = currentLine(iloc+1:Len(currentLine))
          
          read(currentLine, *) ARMdata_struc(n)%stnbd(k)
          
          write(LIS_logunit,*) ARMdata_struc(n)%stn_name(k), &
               ARMdata_struc(n)%stnlat(k), ARMdata_struc(n)%stnlon(k),&
               ARMdata_struc(n)%stnbd(k)
       enddo
    
       call LIS_releaseUnitNumber(ftn)

       ARMdata_struc(n)%udef   = 9999.9

    enddo

    call ESMF_ConfigFindLabel(LIS_config, "ARM objective space attributes file:",rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,objspaceAttribFile(n),rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"ARM number of observation types:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,nobstypes,&
            rc=status)
       call LIS_verify(status)
   
       allocate(obsField(nobstypes))
       allocate(vname(nobstypes))

       ftn=LIS_getNextUnitNumber()
       open(ftn,file=trim(objspaceAttribFile(n)),status='old')
       read(ftn,*)
       
       do i=1,nobstypes
          read(ftn,fmt='(a100)') vname(i)
          obsField(i) = ESMF_FieldCreate(arrayspec=realarrspec, grid=LIS_vecGrid(n), &
               name=trim(vname(i)), rc=status)
          call LIS_verify(status)

          call ESMF_StateAdd(Obj_Space(n),(/obsField(i)/),rc=status)
          call LIS_verify(status)
       enddo
       call LIS_releaseUnitNumber(ftn)
       deallocate(vname)

       call ESMF_TimeIntervalSet(ARMdata_struc(n)%baebbr_ts, s=1800,rc=status)
       call LIS_verify(status, 'Error in timeintervalset: ARM_obsMod ')
       
       call ESMF_TimeIntervalSet(ARMdata_struc(n)%ecor_ts, s=1800,rc=status)
       call LIS_verify(status, 'Error in timeintervalset: ARM_obsMod ')

       call ESMF_TimeIntervalSet(ARMdata_struc(n)%ebbr_ts, s=1800,rc=status)
       call LIS_verify(status, 'Error in timeintervalset: ARM_obsMod ')

       allocate(ARMdata_struc(n)%baebbr_qle_p(ARMdata_struc(n)%n_stns, 48))
       allocate(ARMdata_struc(n)%baebbr_qh_p(ARMdata_struc(n)%n_stns, 48))
       allocate(ARMdata_struc(n)%baebbr_qg_p(ARMdata_struc(n)%n_stns, 48))        
       allocate(ARMdata_struc(n)%baebbr_qle_c(ARMdata_struc(n)%n_stns, 48))
       allocate(ARMdata_struc(n)%baebbr_qh_c(ARMdata_struc(n)%n_stns, 48))
       allocate(ARMdata_struc(n)%baebbr_qg_c(ARMdata_struc(n)%n_stns, 48))        
       allocate(ARMdata_struc(n)%baebbr_tindex(ARMdata_struc(n)%n_stns,48))
       
       allocate(ARMdata_struc(n)%ecor_qle_p(ARMdata_struc(n)%n_stns, 48))
       allocate(ARMdata_struc(n)%ecor_qh_p(ARMdata_struc(n)%n_stns, 48))
       allocate(ARMdata_struc(n)%ecor_qle_c(ARMdata_struc(n)%n_stns, 48))
       allocate(ARMdata_struc(n)%ecor_qh_c(ARMdata_struc(n)%n_stns, 48))
       allocate(ARMdata_struc(n)%ecor_tindex(ARMdata_struc(n)%n_stns,48))

       allocate(ARMdata_struc(n)%ebbr_sfsm_p(ARMdata_struc(n)%n_stns, 48))
       allocate(ARMdata_struc(n)%ebbr_sfst_p(ARMdata_struc(n)%n_stns, 48))
       allocate(ARMdata_struc(n)%ebbr_sfsm_c(ARMdata_struc(n)%n_stns, 48))
       allocate(ARMdata_struc(n)%ebbr_sfst_c(ARMdata_struc(n)%n_stns, 48))

       allocate(ARMdata_struc(n)%ebbr_tindex(ARMdata_struc(n)%n_stns,48))

       ARMdata_struc(n)%ebbr_select = 1
       ARMdata_struc(n)%baebbr_select = 1
       ARMdata_struc(n)%ecor_select = 1
       
       ARMdata_struc(n)%da = -1
       ARMdata_struc(n)%startFlag = .true. 

    enddo
    write(LIS_logunit,*) 'Created the States to hold ARM observations'
    
  end subroutine ARMdata_setup
  
end module ARMdata_module
