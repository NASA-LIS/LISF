!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
module LDT_albedoMod
!BOP
!
! !MODULE: LDT_albedoMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read albedo fraction
!  data. 
!  \subsubsection{Overview}
!  This routines in this module provides routines to read the 
!  albedo fraction climatology data and allows the users to 
!  specify the frequency of climatology (in months). 
!  The climatological data is temporally interpolated  
!  between months to the current simulation date. 
!
! !REVISION HISTORY:
!
!  08 Aug 2005: Sujay Kumar; Initial implementation
!  10 Aug 2012: Kristi Arsenault: Expanded albedo and max snow alb parameters
!  10 Feb 2013: Kristi Arsenault: Modified for CLSM albedo scale factor parameters
!
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use ESMF
  use LDT_coreMod
  use LDT_historyMod
  use LDT_paramDataMod
  use LDT_logMod
  use LDT_paramMaskCheckMod

  implicit none

  PRIVATE

!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LDT_albedo_readParamSpecs
  public :: LDT_albedo_init    !allocates memory for required structures
  public :: LDT_albedo_writeHeader
  public :: LDT_albedo_writeData
!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LDT_albedo_struc

  type, public :: albedo_type_dec

     real             :: alb_gridDesc(20)
     character*50     :: alb_proj
     character*50     :: alb_gridtransform
     character*100    :: albdir
     character*140    :: albfile
     character*20     :: albInterval

     real             :: mxsnoalb_gridDesc(20)
     character*50     :: mxsnoalb_proj
     character*50     :: mxsnoalb_gridtransform
     character*100    :: mxsnoalbfile

   ! Albedo parameters
     type(LDT_paramEntry) :: albedo      ! Climatology-based albedo
     type(LDT_paramEntry) :: mxsnoalb    ! Maximum snow albedo

  end type albedo_type_dec

  type(albedo_type_dec), allocatable :: LDT_albedo_struc(:)


contains

  subroutine LDT_albedo_readParamSpecs
    
    character*100    :: source
    integer          :: rc
    integer          :: n
    
    allocate(LDT_albedo_struc(LDT_rc%nnest))

    call ESMF_ConfigFindLabel(LDT_config,"Albedo data source:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
       call LDT_set_param_attribs(rc,LDT_albedo_struc(n)%albedo,&
            "ALBEDO",source)
    enddo

  ! LSM-required parameter check:
    if( index(LDT_rc%lsm,"Noah") == 1 ) then ! .or. &
!        index(LDT_rc%lsm,"CLSM") == 1 ) then
      if( rc /= 0 ) then
         call LDT_warning(rc,"[WARN] Albedo data source: not defined")
      endif
    elseif( index(LDT_rc%lsm,"CLSM") == 1 ) then
      write(LDT_logunit,*) &
       "[INFO] CLSM F2.5 :: Currently albedo maps are built-in to the model parameter side"
    endif

    call ESMF_ConfigFindLabel(LDT_config,"Max snow albedo data source:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
       call LDT_set_param_attribs(rc,LDT_albedo_struc(n)%mxsnoalb,&
            "MXSNALBEDO",source)
    enddo
  ! LSM-required parameter check:
    if( index(LDT_rc%lsm,"Noah") == 1 .or. &
        index(LDT_rc%lsm,"CLSM") == 1 .or. &
        index(LDT_rc%lsm,"RDHM") == 1 .or. &
        index(LDT_rc%lsm,"SACHTET") == 1 ) then
      if( rc /= 0 ) then
         call LDT_warning(rc,"[WARN] Max snow albedo data source: not defined")
      endif
    endif

  end subroutine LDT_albedo_readParamSpecs

!BOP
! 
! !ROUTINE: LDT_albedo_init
! \label{LDT_albedo_init}
! 
! !INTERFACE:
  subroutine LDT_albedo_init

! !USES:
    use LDT_fileIOMod, only : LDT_readDomainConfigSpecs
    use LDT_logMod,    only : LDT_verify
    use LDT_paramOptCheckMod, only: LDT_albedoOptChecks, &
                       LDT_gridOptChecks

! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! the albedo related datasets.
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[albedosetup](\ref{albedosetup}) \newline
!    calls the registry to invoke the albedo setup methods. 
!  \end{description}
!
!EOP
    implicit none
    integer   :: n
    integer   :: k
    integer   :: rc
    character*3 :: months(12)
    character*3 :: quarters(4)
    data months /'jan','feb','mar','apr','may','jun','jul',&
                 'aug','sep','oct','nov','dec'/
    data quarters /'win','spr','sum','aut'/

    type(LDT_fillopts) :: albedo
    type(LDT_fillopts) :: maxsnalb
    logical            :: alb_select,mxsnoalb_select

    real, allocatable              :: alb_gridDesc(:,:)
    character*50                   :: alb_proj
    character*50,  allocatable     :: alb_gridtransform(:)
    character*100, allocatable     :: albdir(:)
    character*140, allocatable     :: albfile(:)
    character*20,  allocatable     :: albInterval(:)

    real, allocatable              :: mxsnoalb_gridDesc(:,:)
    character*50                   :: mxsnoalb_proj
    character*50,  allocatable     :: mxsnoalb_gridtransform(:)
    character*100, allocatable     :: mxsnoalbfile(:)
! ______________________________________________________________________

    alb_select = .false. 
    do n=1,LDT_rc%nnest
       if(LDT_albedo_struc(n)%albedo%selectOpt.gt.0) then 
          alb_select = .true. 
       endif
    enddo

    mxsnoalb_select = .false. 
    do n=1,LDT_rc%nnest
       if(LDT_albedo_struc(n)%mxsnoalb%selectOpt.gt.0) then 
          mxsnoalb_select = .true. 
       endif
    enddo

    if(alb_select .or. mxsnoalb_select) then
       write(LDT_logunit,*)" - - - - - - - - - - Albedo Parameters - - - - - - - - - - - - -"
    endif

    if ( alb_select ) then
       allocate(alb_gridDesc(LDT_rc%nnest,20))
       allocate(alb_gridtransform(LDT_rc%nnest))         
       allocate(albfile(LDT_rc%nnest))
       allocate(albdir(LDT_rc%nnest))
       allocate(albInterval(LDT_rc%nnest))
    endif

    if ( mxsnoalb_select ) then
       allocate(mxsnoalb_gridDesc(LDT_rc%nnest,20))
       allocate(mxsnoalb_gridtransform(LDT_rc%nnest))         
       allocate(mxsnoalbfile(LDT_rc%nnest))
    endif

    do n = 1,LDT_rc%nnest
       ! Snow-free source and attribs:
       if( alb_select ) then
          call set_albedo_attribs( n, LDT_albedo_struc(n)%albedo%source )
       end if

       ! Max snow-free source and attribs:
       if( mxsnoalb_select ) then
          LDT_albedo_struc(n)%mxsnoalb%vlevels = &
              LDT_albedo_struc(n)%mxsnoalb%num_times

          allocate(LDT_albedo_struc(n)%mxsnoalb%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_albedo_struc(n)%mxsnoalb%vlevels))       
       endif
    enddo   ! end nest loop

! - Read in inputs from config file:

 !- Albedo climatology files:
    if( alb_select ) then

      ! Snow-free albedo interval input:
      call ESMF_ConfigFindLabel(LDT_config,"Albedo climatology interval:",rc=rc)
      do n = 1, LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,albInterval(n),rc=rc)
         call LDT_verify(rc,'Albedo climatology interval: not specified')
         if( trim(albInterval(n)) .ne. "monthly" .and. &
             trim(albInterval(n)) .ne. "quarterly" ) then
           write(LDT_logunit,*) "[ERR] Neither 'monthly' or 'quarterly' albedo interval option specified."
           write(LDT_logunit,*) " Please select one of these two options ... Stopping."
           call LDT_endrun
         endif
         LDT_albedo_struc(n)%albInterval = albInterval(n)

         if( LDT_albedo_struc(n)%albInterval == "monthly" ) then
            LDT_albedo_struc(n)%albedo%num_times = 12
         elseif( LDT_albedo_struc(n)%albInterval == "quarterly" ) then
            LDT_albedo_struc(n)%albedo%num_times = 4
         else
            LDT_albedo_struc(n)%albedo%num_times = 1
         endif

         ! Set climatology interval - dimension:
         LDT_albedo_struc(n)%albedo%vlevels = &
             LDT_albedo_struc(n)%albedo%num_times

         ! Allocate
         allocate(LDT_albedo_struc(n)%albedo%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            LDT_albedo_struc(n)%albedo%vlevels))

         if( albInterval(n) == "monthly" .and. &
             LDT_albedo_struc(n)%albedo%vlevels .ne. 12 .and. &
             LDT_albedo_struc(n)%albedo%selectOpt.eq.1)then
            write(LDT_logunit,*) "[ERR] The 'monthly' albedo interval option should be '12' in the"
            write(LDT_logunit,*) " albedo dimension. Please change to '12' in the code.  Stopping."
            call LDT_endrun
         elseif( albInterval(n) == "quarterly" .and. &
             LDT_albedo_struc(n)%albedo%vlevels .ne. 4 .and.&
             LDT_albedo_struc(n)%albedo%selectOpt.eq.1)then
            write(LDT_logunit,*) "[ERR] The 'quarterly' albedo interval option should be '4' in the"
            write(LDT_logunit,*) " albedo dimension. Please change to '4' in the code.  Stopping."
            call LDT_endrun
         endif
      enddo

      ! Albedo file naming convention, depending on source and interval:
      call ESMF_ConfigFindLabel(LDT_config,"Albedo map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,albdir(n),rc=rc)
         call LDT_verify(rc,'Albedo map: not specified')
         LDT_albedo_struc(n)%albdir = albdir(n)
      enddo

      ! Albedo fill options from ldt.config file:
      albedo%filltype = "none"
      call ESMF_ConfigGetAttribute(LDT_config, albedo%filltype, &
           label="Albedo fill option:",rc=rc)
      call LDT_verify(rc,"Albedo fill option: option not specified in the config file")

      if( albedo%filltype == "neighbor" .or. albedo%filltype == "average" ) then
        call ESMF_ConfigGetAttribute(LDT_config, albedo%fillvalue, &
             label="Albedo fill value:",rc=rc)
        call LDT_verify(rc,"Albedo fill value: option not specified in the config file")

        call ESMF_ConfigGetAttribute(LDT_config, albedo%fillradius, &
             label="Albedo fill radius:",rc=rc)
        call LDT_verify(rc,"Albedo fill radius: option not specified in the config file")
      elseif( albedo%filltype == "none" ) then
         write(LDT_logunit,*) "[INFO] 'NONE' Parameter-Mask Agreement Option Selected for Albedo"
      else
         write(LDT_logunit,*) "[ERR] Fill option for Albedo is not valid: ",trim(albedo%filltype)
         write(LDT_logunit,*) " Please select one of these:  none, neighbor or average "
         write(LDT_logunit,*) " Programming stopping ..."
         call LDT_endrun
      end if

    ! Read in Albedo map options:
      call ESMF_ConfigGetAttribute(LDT_config,alb_proj,&
           label="Albedo map projection:",rc=rc)
      call LDT_verify(rc,'Albedo projection: option not specified in the config file')
      LDT_albedo_struc(:)%alb_proj = alb_proj

      call ESMF_ConfigFindLabel(LDT_config,"Albedo spatial transform:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,alb_gridtransform(n),&
              rc=rc)
         call LDT_verify(rc,'Albedo spatial transform: option not specified in the config file')
         LDT_albedo_struc(n)%alb_gridtransform = alb_gridtransform(n)
      enddo

    ! Don't need to read in "Native" grid extents/resolution, just for LIS inputs
      do n=1,LDT_rc%nnest
         if( index(LDT_albedo_struc(n)%albedo%source,"Native").eq. 0  .and. &
             index(LDT_albedo_struc(n)%albedo%source,"CONSTANT").eq. 0  ) then
            call LDT_readDomainConfigSpecs("Albedo", alb_proj, alb_gridDesc)
            if( alb_proj == "latlon" ) then
              call LDT_gridOptChecks( n, "Albedo", &
                   alb_gridtransform(n), &
                   alb_proj, alb_gridDesc(n,9) )
              LDT_albedo_struc(n)%alb_gridDesc = alb_gridDesc(n,:)
           endif
        endif
      end do

    ! Set units and full names:
      do n=1,LDT_rc%nnest
         LDT_albedo_struc(n)%albedo%units="-"
         call setAlbedoParmsFullnames( n, "albedo", &
                 LDT_albedo_struc(n)%albedo%source )
      enddo

    end if !- End albedo clim read

    
  ! Max snow albedo:
    if(mxsnoalb_select) then

       call ESMF_ConfigFindLabel(LDT_config,"Max snow albedo map:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,mxsnoalbfile(n),rc=rc)
          call LDT_verify(rc,'Max snow albedo map: not specified')
          LDT_albedo_struc(n)%mxsnoalbfile = mxsnoalbfile(n)
       enddo
       call ESMF_ConfigGetAttribute(LDT_config,mxsnoalb_proj,&
            label="Max snow albedo map projection:",rc=rc)
       call LDT_verify(rc,'Max snow albedo projection: option not specified in the config file')
       LDT_albedo_struc(:)%mxsnoalb_proj = mxsnoalb_proj

       call ESMF_ConfigFindLabel(LDT_config,"Max snow albedo spatial transform:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,mxsnoalb_gridtransform(n),&
               rc=rc)
          call LDT_verify(rc,'Max snow albedo spatial transform: option not specified in the config file')
          LDT_albedo_struc(n)%mxsnoalb_gridtransform = mxsnoalb_gridtransform(n)
       enddo

     ! Don't need to read in "Native" grid extents/resolution, just for LIS inputs
       do n=1,LDT_rc%nnest
          if( index(LDT_albedo_struc(n)%mxsnoalb%source,"Native").eq.0  .and. &
              index(LDT_albedo_struc(n)%mxsnoalb%source,"CONSTANT").eq.0) then
            call LDT_readDomainConfigSpecs("Max snow albedo", &
                 mxsnoalb_proj, &
                 mxsnoalb_gridDesc)
            LDT_albedo_struc(n)%mxsnoalb_gridDesc = mxsnoalb_gridDesc(n,:)
            if( LDT_albedo_struc(n)%mxsnoalb_proj == "latlon" ) then
               call LDT_gridOptChecks(n,"Max snow albedo", &
                    mxsnoalb_gridtransform(n), &
                    mxsnoalb_proj, mxsnoalb_gridDesc(n,9))
               LDT_albedo_struc(n)%mxsnoalb_gridDesc = mxsnoalb_gridDesc(n,:)
            endif
          endif
       enddo

    ! Set units and full names:
      do n=1,LDT_rc%nnest
         LDT_albedo_struc(n)%mxsnoalb%units="-"
         call setAlbedoParmsFullnames( n, "maxsnowalb", &
                 LDT_albedo_struc(n)%mxsnoalb%source )
      enddo

       maxsnalb%filltype = "none"
       call ESMF_ConfigGetAttribute(LDT_config, maxsnalb%filltype, &
            label="Max snow albedo fill option:",rc=rc)
       call LDT_verify(rc,"Max snow albedo fill option: option not specified in the config file")

       if( maxsnalb%filltype == "neighbor" .or. maxsnalb%filltype == "average" ) then
         call ESMF_ConfigGetAttribute(LDT_config, maxsnalb%fillvalue, &
              label="Max snow albedo fill value:",rc=rc)
         call LDT_verify(rc,"Max snow albedo fill value: option not specified in the config file")
 
         call ESMF_ConfigGetAttribute(LDT_config, maxsnalb%fillradius, &
              label="Max snow albedo fill radius:",rc=rc)
         call LDT_verify(rc,"Max snow albedo fill radius: option not specified in the config file")
       elseif( maxsnalb%filltype == "none" ) then
         write(LDT_logunit,*) "[INFO] Parameter-Mask Agreement Option Selected for Max snow albedo"
       else
         write(LDT_logunit,*) "[ERR] Fill option for Max snow albedo is not valid: ",trim(albedo%filltype)
         write(LDT_logunit,*) "  Please select one of these:  none, neighbor or average "
         write(LDT_logunit,*) "  Programming stopping ..."
         call LDT_endrun
       end if

    endif    ! Max snow albedo selected

!-- Read in albedo files: 
    do n = 1, LDT_rc%nnest

    !- Maximum snow albedo parameter:
       if(LDT_albedo_struc(n)%mxsnoalb%selectOpt.eq.1) then

          call LDT_albedoOptChecks( "Max snow albedo", mxsnoalb_proj, &
               mxsnoalb_gridtransform(n) )

          write(LDT_logunit,*) 'Reading '//trim(mxsnoalbfile(n))
          call readmxsnoalb( trim(LDT_albedo_struc(n)%mxsnoalb%source)//char(0),&
               n,LDT_albedo_struc(n)%mxsnoalb%value(:,:,1) )
          write(LDT_logunit,*) 'Done reading '//trim(mxsnoalbfile(n))

        ! Fill where parameter values are missing compared to land/water mask:
          if( maxsnalb%filltype == "neighbor" .or. &
               maxsnalb%filltype == "average" ) then
             write(LDT_logunit,*) "Checking/filling mask values for: ", &
                  trim(LDT_albedo_struc(n)%mxsnoalb%short_name)
             write(fill_logunit,*) "Checking/filling mask values for: ", &
                  trim(LDT_albedo_struc(n)%mxsnoalb%short_name)
             maxsnalb%watervalue = LDT_rc%udef  ! LIS-based NCEP params
             
             call LDT_contIndivParam_Fill( n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
                  mxsnoalb_gridtransform(n),                                &
                  LDT_albedo_struc(n)%mxsnoalb%num_bins,                 &
                  LDT_albedo_struc(n)%mxsnoalb%value, maxsnalb%watervalue, &
                  LDT_LSMparam_struc(n)%landmask2%value,                    &
                  maxsnalb%filltype, maxsnalb%fillvalue, maxsnalb%fillradius )
          endif
       endif
       
    !- Snow-free albedo parameter:
       if(LDT_albedo_struc(n)%albedo%selectOpt.eq.1) then

          call LDT_albedoOptChecks( "Albedo", alb_proj, &
               alb_gridtransform(n) )

          if( albInterval(n) == "monthly" ) then     
             LDT_rc%monthlyData(n) = .true.
          elseif( albInterval(n) == "quarterly") then 
             LDT_rc%quarterlyData(n) = .true.
          endif

       !- Read in files:
          do k = 1, LDT_albedo_struc(n)%albedo%vlevels

             if(albInterval(n).eq."monthly") then      
               if(trim(LDT_albedo_struc(n)%albedo%source).eq."NCEP_LIS") then
                  albfile(n) = trim(albdir(n))//'.'//&
                                      trim(months(k))//'.1gd4r'
               elseif(trim(LDT_albedo_struc(n)%albedo%source).eq."NCEP_Native") then
                  albfile(n) = trim(albdir(n))//'_'//&
                                      trim(months(k))//'.asc'
               elseif(trim(LDT_albedo_struc(n)%albedo%source).eq."CONSTANT") then
                 albfile(n) = trim(albdir(n))
               endif
             elseif(albInterval(n).eq."quarterly") then
               if(trim(LDT_albedo_struc(n)%albedo%source).eq."NCEP_LIS") then
                 albfile(n) = trim(albdir(n))//'.'//&
                                     trim(quarters(k))//'.1gd4r'
               elseif(trim(LDT_albedo_struc(n)%albedo%source).eq."NCEP_Native") then
                 albfile(n) = trim(albdir(n))
               elseif(trim(LDT_albedo_struc(n)%albedo%source).eq."CONSTANT") then
                 albfile(n) = trim(albdir(n))
               endif
             endif
             LDT_albedo_struc(n)%albfile = albfile(n)
                            
             write(LDT_logunit,*) "Reading "//trim(albfile(n))

          !- Account for original NCEP-Native quarterly (ASCII) files separately: 
             if( LDT_albedo_struc(n)%albedo%source.eq."NCEP_Native" .and. &
                 albInterval(n).eq."quarterly") then 
                call readalbedo( "NCEP_NativeQtr"//char(0),&
                                 n,LDT_albedo_struc(n)%albedo%value(:,:,k))
                write(LDT_logunit,*) "Done reading "//trim(albfile(n))
             !- Exit time-dimension loop here, if quarterly and NCEP_Native file read in: 
                exit

             else  ! All other albedo files:
                call readalbedo(trim(LDT_albedo_struc(n)%albedo%source)//char(0),&
                                n,LDT_albedo_struc(n)%albedo%value(:,:,k))
             endif
             write(LDT_logunit,*) "Done reading "//trim(albfile(n))
          enddo

        ! Fill where parameter values are missing compared to land/water mask:
          if( albedo%filltype == "neighbor" .or. &
              albedo%filltype == "average" ) then
            write(LDT_logunit,*) "Checking/filling mask values for: ", &
                              trim(LDT_albedo_struc(n)%albedo%short_name)
            write(fill_logunit,*) "Checking/filling mask values for: ", &
                              trim(LDT_albedo_struc(n)%albedo%short_name)
            albedo%watervalue = LDT_rc%udef
            call LDT_contIndivParam_Fill( n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
                  alb_gridtransform(n),                          &
                  LDT_albedo_struc(n)%albedo%vlevels,              &
                  LDT_albedo_struc(n)%albedo%value, albedo%watervalue, &
                  LDT_LSMparam_struc(n)%landmask2%value,              &
                  albedo%filltype, albedo%fillvalue, albedo%fillradius )
          endif
       endif
    enddo

    if ( alb_select ) then
       deallocate(alb_gridDesc)
       deallocate(alb_gridtransform)         
       deallocate(albfile)
       deallocate(albdir)
       deallocate(albInterval)
    endif

    if ( mxsnoalb_select ) then
       deallocate(mxsnoalb_gridDesc)
       deallocate(mxsnoalb_gridtransform)         
       deallocate(mxsnoalbfile)
    endif


  end subroutine LDT_albedo_init


  subroutine LDT_albedo_writeHeader(n,ftn,dimID,monthID,qid)

    integer    :: n 
    integer    :: ftn
    integer    :: dimID(3)
    integer    :: monthID
    integer    :: qID

    integer    :: t_dimID(3)

    t_dimID(1) = dimID(1)
    t_dimID(2) = dimID(2)

 !- Snow-free albedo:
    if( LDT_albedo_struc(n)%albedo%selectOpt.eq.1 ) then
       if( LDT_albedo_struc(n)%albInterval.eq."monthly" ) &   ! monthly
          t_dimID(3) = monthID  
       if( LDT_albedo_struc(n)%albInterval.eq."quarterly" ) &   ! quarterly
          t_dimID(3) = qID  

       call LDT_writeNETCDFdataHeader(n,ftn,t_dimID,&
            LDT_albedo_struc(n)%albedo)

       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"ALBEDO_DATA_INTERVAL", &
            LDT_albedo_struc(n)%albInterval))
    endif

 !- Max snow albedo:
    if( LDT_albedo_struc(n)%mxsnoalb%selectOpt.eq.1 ) &
       call LDT_writeNETCDFdataHeader(n,ftn,t_dimID,&
            LDT_albedo_struc(n)%mxsnoalb)


  end subroutine LDT_albedo_writeHeader


  subroutine LDT_albedo_writeData(n,ftn)

    integer    :: n 
    integer    :: ftn

    if( LDT_albedo_struc(n)%albedo%selectOpt.eq.1 ) &
      call LDT_writeNETCDFdata(n,ftn,LDT_albedo_struc(n)%albedo)

    if( LDT_albedo_struc(n)%mxsnoalb%selectOpt.eq.1 ) &
      call LDT_writeNETCDFdata(n,ftn,LDT_albedo_struc(n)%mxsnoalb)


  end subroutine LDT_albedo_writeData

end module LDT_albedoMod
