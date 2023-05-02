!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_DAobsDataMod
!BOP
!
!  !MODULE: LDT_DAobsDataMod
! 
!  !DESCRIPTION: 
!   This module is used to define the metadata associated with 
!   observational data. The list of supported variables are defined
!   based ont the ALMA specification. 
!
!   \textsl{http://www.lmd.jussieu.fr/ALMA/}
!   
!  !REVISION HISTORY: 
!  2 Oct 2008    Sujay Kumar  Initial Specification
!  2 Dec 2021:   Mahdi Navari; modified to compute CDF for precipitation
! !USES: 

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS: 
!-----------------------------------------------------------------------------
  public :: LDT_DAobsEntryInit
  public :: LDT_initializeDAobsEntry
  public :: LDT_tavgDAobsData
  public :: LDT_resetDAobsData
  public :: LDT_logSingleDAobs

!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
! MOC - Model Output Convention
!-----------------------------------------------------------------------------
  public :: LDT_DAobsData
  public :: LDT_DAobsDataPtr
  public :: LDT_DA_MOC_SWE       
  public :: LDT_DA_MOC_SNOWDEPTH
  public :: LDT_DA_MOC_SOILMOIST 
  public :: LDT_DA_MOC_SOILTEFF   !Y.Kwon
  public :: LDT_DA_MOC_TWS
  public :: LDT_DA_MOC_VOD
  public :: LDT_DA_MOC_LAI
  public :: LDT_DA_MOC_GVF   !Y.Kwon
  public :: LDT_DA_MOC_COUNT
  public :: LDT_DA_MOC_TOTALPRECIP 
!  public :: LDT_MOC_GRIB_COUNT

   ! ALMA ENERGY BALANCE COMPONENTS
  integer, parameter :: LDT_DA_MOC_SWE        = 1
  integer, parameter :: LDT_DA_MOC_SNOWDEPTH  = 2
  integer, parameter :: LDT_DA_MOC_SOILMOIST  = 3
  integer, parameter :: LDT_DA_MOC_TWS        = 4
  integer, parameter :: LDT_DA_MOC_VOD        = 5
  integer, parameter :: LDT_DA_MOC_LAI        = 6
  integer, parameter :: LDT_DA_MOC_GVF        = 7    !Y.Kwon
  integer, parameter :: LDT_DA_MOC_SOILTEFF   = 8    !Y.Kwon
  integer, parameter :: LDT_DA_MOC_TOTALPRECIP= 9
  ! READ ABOVE NOTE ABOUT SPECIAL CASE INDICES
  integer, parameter :: LDT_DA_MOC_COUNT      = 9     !Y.Kwon
  ! Add the special cases.  LDT_MOC_GRIB_COUNT should be used only in
   ! LDT_gribMod.F90.
!  integer, parameter :: LDT_MOC_GRIB_COUNT = 100
  
!  real, parameter :: LDT_MOC_MAX_NUM =  999999.0
!  real, parameter :: LDT_MOC_MIN_NUM = -999999.0

!EOP
  type, public :: LDT_DAmetadataEntry
     character*100         :: long_name
     character*20          :: standard_name
     character*20          :: units
     integer               :: nunits
     character*1           :: format          ! (scientific - E, else - F)
     character*20          :: dir
     real                  :: valid_min
     real                  :: valid_max
     integer               :: vlevels         ! Number of classes or time points
     integer               :: num_bins        ! Number of bins 
     integer               :: gribId
     integer               :: gribSF          ! grib scale factor
     integer               :: gribCat
     integer               :: timeAvgOpt
     integer               :: selectOpt
     integer               :: selectStats      !whether to output stats
     integer               :: minMaxOpt
     integer               :: stdOpt
     integer               :: varID(6) !varID(5)  !YK
     character*20, allocatable :: unittypes(:)
     integer, allocatable      :: count(:,:)
     real, allocatable         :: value(:,:) 
  end type LDT_DAmetadataEntry

  type, public :: output_meta

     type(LDT_DAmetadataEntry) :: swe          ! Snow water equivalent (kg/m2)
     type(LDT_DAmetadataEntry) :: snowdepth
     type(LDT_DAmetadataEntry) :: soilmoist
     type(LDT_DAmetadataEntry) :: teff         ! effective soil temperature (Y.Kwon)
     type(LDT_DAmetadataEntry) :: vod
     type(LDT_DAmetadataEntry) :: lai
     type(LDT_DAmetadataEntry) :: totalprecip
     type(LDT_DAmetadataEntry) :: gvf          !Y.Kwon

  end type output_meta

  type, public :: dep
     type(LDT_DAmetadataEntry), pointer :: dataEntryPtr
  end type dep

  type, public ::  obs_list_dec
     type(LDT_DAmetadataEntry) :: swe_obs
     type(LDT_DAmetadataEntry) :: snowdepth_obs
     type(LDT_DAmetadataEntry) :: soilmoist_obs
     type(LDT_DAmetadataEntry) :: teff_obs      !Y.Kwon
     type(LDT_DAmetadataEntry) :: tws_obs
     type(LDT_DAmetadataEntry) :: vod_obs
     type(LDT_DAmetadataEntry) :: lai_obs
     type(LDT_DAmetadataEntry) :: totalprecip_obs
     type(LDT_DAmetadataEntry) :: gvf_obs      !Y.Kwon
  end type obs_list_dec

  type, public :: obsdep
     type(LDT_DAmetaDataEntry), pointer :: dataEntryPtr
  end type obsdep

  type(obs_list_dec),     allocatable :: LDT_DAobsData(:)
  type(obsdep),           pointer :: LDT_DAobsDataPtr(:,:)
contains


!BOP
! 
! !ROUTINE: allocate_dataEntry
! \label{allocate_dataEntry}
! 
! !INTERFACE: 
  subroutine allocate_dataEntry(dataEntry, nunits,nsize,unittypes)
    
    implicit none
! !ARGUMENTS: 
    type(LDT_DAmetadataEntry) :: dataEntry
    integer                 :: nunits
    integer                 :: nsize
    character(len=*)        :: unittypes(nunits)
! 
! !DESCRIPTION: 
!  This routine initializes the datastructures required for the 
!  specified output variable. 
!EOP
    integer                 :: i 

    if(dataEntry%selectOpt.ne.0) then 
       allocate(dataEntry%value(nsize,dataEntry%vlevels))
       allocate(dataEntry%count(nsize,dataEntry%vlevels))

       dataEntry%value = 0 
       dataEntry%count = 0 

       dataEntry%nunits = nunits
       allocate(dataEntry%unittypes(nunits))
       
       do i=1,nunits
          dataEntry%unittypes(i) = trim(unittypes(i))
       enddo
    endif
  end subroutine allocate_dataEntry


!BOP
! !ROUTINE: LDT_DAobsEntryInit
! \label{LDT_DAobsEntryInit}
! 
! !INTERFACE: 
  subroutine LDT_DAobsEntryInit(i, nsize)
! !ARGUMENTS: 
    integer,  intent(IN)   :: i 
    integer,  intent(IN)   :: nsize
! 
! !DESCRIPTION: 
!  This routine initializes the required arrays to hold the selected list of 
!   LSM variables
!
!   The arguments are: 
!   \begin{description}
!    \item[i]  index of the observational plugin \newline
!    \item[nsize]  size of the tilespace \newline
!   \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!    \item[register\_obsDataEntry] (\ref{register_obsDataEntry}) \newline
!      registers the data structures related to the variable being output
!   \end{description}
!EOP
    call register_obsDataEntry(i,LDT_DA_MOC_SWE ,&
         LDT_DAobsData(i)%swe_obs,1,nsize,(/"kg/m2"/))
    call register_obsDataEntry(i,LDT_DA_MOC_SNOWDEPTH ,&
         LDT_DAobsData(i)%snowdepth_obs,1,nsize,(/"m"/))
    call register_obsDataEntry(i,LDT_DA_MOC_SOILMOIST ,&
         LDT_DAobsData(i)%soilmoist_obs,2,nsize,(/"kg/m2", "m3/m3"/))
    call register_obsDataEntry(i,LDT_DA_MOC_SOILTEFF ,&
         LDT_DAobsData(i)%teff_obs,1,nsize,(/"K"/))     !Y.Kwon
    call register_obsDataEntry(i,LDT_DA_MOC_TWS ,&
         LDT_DAobsData(i)%tws_obs,1,nsize,(/"mm"/))
    call register_obsDataEntry(i,LDT_DA_MOC_VOD ,&
         LDT_DAobsData(i)%vod_obs,1,nsize,(/"-"/))
    call register_obsDataEntry(i,LDT_DA_MOC_LAI ,&
         LDT_DAobsData(i)%lai_obs,1,nsize,(/"-"/))
    call register_obsDataEntry(i,LDT_DA_MOC_TOTALPRECIP ,&
         LDT_DAobsData(i)%totalprecip_obs,1,nsize,(/"kg/m2"/))
    call register_obsDataEntry(i,LDT_DA_MOC_GVF ,&
         LDT_DAobsData(i)%gvf_obs,1,nsize,(/"-"/))    !Y.Kwon
  end subroutine LDT_DAobsEntryInit

!BOP
! !ROUTINE: LDT_initializeDAobsEntry
! \label{LDT_initializeDAobsEntry}
! 
! !INTERFACE: 
  subroutine LDT_initializeDAobsEntry(dataEntry, unit, timeAvgOpt, vlevel)

    implicit none
! !ARGUMENTS: 
    type(LDT_DAmetadataEntry)      :: dataEntry
    character(len=*)             :: unit
    integer                      :: timeAvgOpt
    integer                      :: vlevel 
!
! !DESCRIPTION: 
!    This subroutine initializes the obs data structure with the 
!    values of units, vertical levels and options for selection and 
!    time averaging.  
!  
!EOP
    dataEntry%selectOpt  = 1
    dataEntry%units      = trim(unit)
    dataEntry%timeAvgOpt = timeAvgOpt
    dataEntry%vlevels    = vlevel 

  end subroutine LDT_initializeDAobsEntry
  
!BOP
! 
! !ROUTINE: register_obsDataEntry
! \label{register_obsDataEntry}
! 
! !INTERFACE: 
  subroutine register_obsDataEntry(i,ldt_moc_index,dataEntry, &
       nunits,nsize,unittypes)
    
    implicit none
! !ARGUMENTS: 
    integer                 :: i 
    integer                 :: ldt_moc_index
    type(LDT_DAmetadataEntry), target :: dataEntry
    integer                 :: nunits
    integer                 :: nsize
    character(len=*)        :: unittypes(nunits)
! 
! !DESCRIPTION: 
!  This routine initializes the datastructures required for the 
!  specified observational output variable. 
!EOP

    LDT_DAobsDataPtr(i,ldt_moc_index)%dataEntryPtr => dataEntry
    call allocate_obsDataEntry(dataEntry, nunits, nsize, unittypes)

  end subroutine register_obsDataEntry

  
!BOP
! 
! !ROUTINE: allocate_obsDataEntry
! \label{allocate_obsDataEntry}
! 
! !INTERFACE: 
  subroutine allocate_obsDataEntry(dataEntry, nunits,nsize,unittypes)
    
    implicit none
! !ARGUMENTS: 
    type(LDT_DAmetadataEntry) :: dataEntry
    integer                 :: nunits
    integer                 :: nsize
    character(len=*)        :: unittypes(nunits)
! 
! !DESCRIPTION: 
!  This routine initializes the datastructures required for the 
!  specified observational output variable. 
!EOP
    integer                 :: i 

    if(dataEntry%selectOpt.ne.0) then 
       allocate(dataEntry%value(nsize,dataEntry%vlevels))
       allocate(dataEntry%count(nsize,dataEntry%vlevels))
       dataEntry%nunits = nunits
       allocate(dataEntry%unittypes(nunits))
       
       dataEntry%value = 0 
       dataEntry%count = 0 

       do i=1,nunits
          dataEntry%unittypes(i) = trim(unittypes(i))
       enddo
    endif
  end subroutine allocate_obsDataEntry

!BOP
! 
! !ROUTINE: LDT_tavgDAObsData
! \label{LDT_tavgDAObsData}
! 
! !INTERFACE: 
  subroutine LDT_tavgDAObsData(n)
! !USES: 
    use LDT_coreMod, only : LDT_rc

! !DESCRIPTION: 
!   This routine invokes the calls to compute temporal averages of 
!   desired set of observational output variables, based on the specified 
!   temporal averaging frequency
!  
!   The routines invoked are: 
!   \begin{description}
!    \item[tavgSingleObs](\ref{tavgSingleObs})
!     computes the temporal average for a single observational variable
!   \end{description}
!EOP
    implicit none
    integer,     intent(in) :: n 
    integer  :: rc
    integer  :: i, index

    if(mod(float(LDT_rc%hr)*3600+60*float(LDT_rc%mn)+float(LDT_rc%ss),&
            LDT_rc%tavgInterval).eq.0) then   
       LDT_rc%computeFlag = .true.
    else
       LDT_rc%computeFlag = .false. 
    endif

    if(LDT_rc%computeFlag) then 
       do index=1, LDT_DA_MOC_COUNT             
          call tavgSingleObs(n,LDT_DAobsDataPtr(n,index)%dataEntryPtr)
       enddo
    endif

  end subroutine LDT_tavgDAObsData

!BOP
!
! !ROUTINE: tavgSingleObs
! \label{tavgSingleObs}
! 
! !INTERFACE:
  subroutine tavgSingleObs(n,dataEntry)
! !USES: 
    use LDT_coreMod, only      : LDT_rc, LDT_domain
    use LDT_logMod,  only      : LDT_logunit, LDT_endrun

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n 
    integer :: ftn
    integer :: form
    type(LDT_DAmetadataEntry) :: dataEntry
! 
! !DESCRIPTION: 
!   This subroutine temporally averages the accumulated observational 
!   data for a single variable. 
!EOP
    integer :: unit_id
    integer :: k,i,c,r,gid

    if(dataEntry%selectOpt.eq.1) then 
       do k=1,dataEntry%vlevels
          do r=1,LDT_rc%lnr(n)
             do c=1,LDT_rc%lnc(n)
                if(LDT_domain(n)%gindex(c,r).ne.-1) then 
                   gid = LDT_domain(n)%gindex(c,r)
                   if(dataEntry%count(gid,k).ne.0) then 
                      dataEntry%value(gid,k) = &
                           dataEntry%value(gid,k)/dataEntry%count(gid,k)
                   endif
                endif
             enddo
          enddo
       enddo
    endif
  end subroutine tavgSingleObs

!BOP
! 
! !ROUTINE: LDT_resetDAobsData
! \label{LDT_resetDAobsData}
! 
! !INTERFACE: 
  subroutine LDT_resetDAobsData(n)
! !USES: 
    use LDT_coreMod, only : LDT_rc
!
! !DESCRIPTION: 
!   This routine reinitializes the data structures that hold the observational
!   data
! 
!   The routines invoked are: 
!   \begin{description}
!    \item[resetSingleObs](\ref{resetSingleObs})
!     resets the datastructures for a single variable
!   \end{description}
!EOP
    implicit none
    integer, intent(in) :: n 
    integer  :: rc
    integer  :: i, index

    do index=1, LDT_DA_MOC_COUNT             
       call resetSingleObs(n, LDT_DAobsDataPtr(n,index)%dataEntryPtr)
    enddo

  end subroutine LDT_resetDAobsData
  
!BOP
! !ROUTINE: resetSingleObs
! \label{resetSingleObs}
! 
! !INTERFACE: 
  subroutine resetSingleObs(n, dataEntry)

    implicit none 
! !ARGUMENTS: 
    type(LDT_DAmetadataEntry) :: dataEntry
! 
! !DESCRIPTION: 
!  This routine resets the data structures that hold the observational 
!  data and the temporal averaging counters
!EOP
    integer, intent(in) :: n 
    integer                 :: k 

    if(dataEntry%selectOpt.eq.1) then 
       do k=1,dataEntry%vlevels
          dataEntry%value(:,k) = 0 
          dataEntry%count(:,k) = 0 
       enddo      
    endif
  end subroutine resetSingleObs

!BOP
!
! !ROUTINE: LDT_logSingleDAobs
! \label{LDT_logSingleDAobs}
! 
! !INTERFACE:
  subroutine LDT_logSingleDAobs(n, dataEntry, value,vlevel)
! !USES: 
    use LDT_coreMod, only     : LDT_rc, LDT_domain
    use LDT_logMod,  only     : LDT_logunit, LDT_endrun

    implicit none
! !ARGUMENTS: 
    integer, intent(in)        :: n 
    type(LDT_DAmetadataEntry)  :: dataEntry
    real                       :: value(LDT_rc%lnc(n), LDT_rc%lnr(n))
    integer, optional          :: vlevel
! 
! !DESCRIPTION: 
!  This subroutine maps the processed observations onto the LDT data
!  structures for future temporal averaging and comparisons. The 
!  data are also filtered using the specified external mask. 
!
!EOP
    integer :: form
    integer :: k,i,c,r,gid

    if(.not.present(vlevel)) then 
       k = 1
    else
       k = vlevel
    endif

    if(dataEntry%selectOpt.eq.1) then 
       do r=1,LDT_rc%lnr(n)
          do c=1,LDT_rc%lnc(n)
!             if(trim(dataEntry%standard_name).eq."SoilMoist") then                
!                if(c.eq.61.and.r.eq.91) print*, value(c,r), dataEntry%value(c,r,k)
!             endif
             if(LDT_domain(n)%datamask(c,r).eq.1) then 
                if(value(c,r).ne.LDT_rc%udef) then
                   if(LDT_domain(n)%gindex(c,r).ne.-1) then 
                      gid = LDT_domain(n)%gindex(c,r)
                      dataEntry%value(gid,k) = dataEntry%value(gid,k) + value(c,r)
                      dataEntry%count(gid,k) = dataEntry%count(gid,k) + 1
                   endif
                endif
             endif
          enddo
       enddo
    endif
  end subroutine LDT_logSingleDAobs
end module LDT_DAobsDataMod
