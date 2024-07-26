!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_climateParmsMod
!BOP
!
! !MODULE: LDT_climateParmsMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read forcing 
!   climatology data (e.g., can be used for forcing downscaling, etc.). 
!
!  \subsubsection{Overview}
!  The routines in this module provide capabilities to read the 
!  climatology data and allows the users to 
!  specify the frequency of climatology (in months). 
!  The climatological data is temporally interpolated  
!  between months to the current simulation date. 
!
! !REVISION HISTORY:
!
!  08 Aug 2005: Sujay Kumar; Initial implementation
!  08 Oct 2012: Kristi Arsenault; Expanded for forcing climatology datasets
!  28 Jun 2022: Eric Kemp; Added NAFPA GFS and GALWEM background datasets.
!
  use ESMF
  use LDT_coreMod
  use LDT_historyMod
  use LDT_metforcingParmsMod
  use LDT_paramDataMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_paramMaskCheckMod

  implicit none

  PRIVATE

!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LDT_climate_readParamSpecs
  public :: LDT_climateparms_init    !allocates memory for required structures
  public :: LDT_climateparms_writeHeader
  public :: LDT_climateparms_writeData

!------------------------------------------------------------------------------
! !PUBLIC TYPE:
!------------------------------------------------------------------------------
  public :: LDT_climate_struc

  type, public :: climate_type_dec
     character*50  :: clim_gridtransform
     character*50  :: clim_gridtransform2

     character(len=LDT_CONST_PATH_LEN) :: climpptdir
     character(len=LDT_CONST_PATH_LEN) :: climpptdir2
     character(len=LDT_CONST_PATH_LEN) :: climtmindir
     character(len=LDT_CONST_PATH_LEN) :: climtmaxdir
     integer       :: output_climppt_ratio

     character*20  :: climpptInterval
     character*20  :: climtminInterval
     character*20  :: climtmaxInterval

     character(len=LDT_CONST_PATH_LEN) :: climpptfile
     character(len=LDT_CONST_PATH_LEN) :: climpptfile2
     character(len=LDT_CONST_PATH_LEN) :: climtmaxfile
     character(len=LDT_CONST_PATH_LEN) :: climtminfile
     character(len=LDT_CONST_PATH_LEN) :: climelevfile

     integer :: climpptimonth

     ! Climate parameters
     type(LDT_paramEntry) :: climppt     ! Precipitation climatology (LIS-domain)
     type(LDT_paramEntry) :: climppt2     ! Precipitation climatology (LIS-domain)
     type(LDT_paramEntry) :: climtmin    ! Min. temp. climatology
     type(LDT_paramEntry) :: climtmax    ! Min. temp. climatology
     type(LDT_paramEntry) :: climelev    ! Forcing climatology elev (LIS-domain)

  end type climate_type_dec

  type(climate_type_dec), allocatable :: LDT_climate_struc(:)

contains

  subroutine LDT_climate_readParamSpecs
    
    character*100   :: source
    integer         :: rc,output_ratio
    integer         :: n

    allocate(LDT_climate_struc(LDT_rc%nnest))
    output_ratio = 0
    call ESMF_ConfigFindLabel(LDT_config,"PPT climatology data source 2:",rc=rc)
    ! EMK...Check source of second climatology.  Don't just assume NLDAS
    if (rc == 0) then
       output_ratio =1
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
          call LDT_set_param_attribs(rc,LDT_climate_struc(n)%climppt2,&
               "PPT_ratio", source)
       enddo
    end if
    call ESMF_ConfigFindLabel(LDT_config,"PPT climatology data source:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
       if (output_ratio > 0) then
         LDT_climate_struc(n)%output_climppt_ratio = 1
         call LDT_set_param_attribs(rc,LDT_climate_struc(n)%climppt,&
            "PPT_ratio",source)
       else
         call LDT_set_param_attribs(rc,LDT_climate_struc(n)%climppt,&
            "PPTCLIM",source)
       end if
    enddo

    if( LDT_rc%nmetforc > 0 .and. rc /= 0 ) then
       call LDT_warning(rc,"PPT climatology data source: not defined")
    endif

! Climate Tmin / Tmax included here - but currently not fully supported.
! Tmin:
    call ESMF_ConfigFindLabel(LDT_config,"Tmin climatology data source:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
       call LDT_set_param_attribs(rc,LDT_climate_struc(n)%climtmin,&
            "CLIMTMIN",source)
    enddo
    if( LDT_rc%nmetforc > 0 .and. rc /= 0 ) then
       call LDT_warning(rc,"Tmin climatology data source: not defined")
    endif
! Tmax:
    call ESMF_ConfigFindLabel(LDT_config,"Tmax climatology data source:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
       call LDT_set_param_attribs(rc,LDT_climate_struc(n)%climtmax,&
            "CLIMTMAX",source)
    enddo
    if( LDT_rc%nmetforc > 0 .and. rc /= 0 ) then
       call LDT_warning(rc,"Tmax climatology data source: not defined")
    endif

  end subroutine LDT_climate_readParamSpecs

!BOP
! 
! !ROUTINE: LDT_climateparms_init
! \label{LDT_climateparms_init}
! 
! !INTERFACE:
  subroutine LDT_climateparms_init

! !USES:
    use LDT_fileIOMod
    use LDT_logMod,    only : LDT_verify
    use LDT_paramOptCheckMod, only: LDT_climateOptChecks, &
                       LDT_gridOptChecks

! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! forcing climatology datasets. 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[climpptsetup](\ref{climpptsetup}) \newline
!    calls the registry to invoke the climppt setup methods. 
!  \end{description}
!
!EOP
    implicit none
    integer   :: n, k, m, ic, ir
    integer   :: rc
    character*3 :: months(12), mon2d(12)
    data months /'jan','feb','mar','apr','may','jun','jul','aug',&
                 'sep','oct','nov','dec'/
    data mon2d /'1','2','3','4','5','6','7','8','9','10','11','12'/

    type(LDT_fillopts) :: climppt
    logical            :: climppt_select
    logical            :: climtmin_select
    logical            :: climtmax_select
    real :: ratio

    external :: setClimateParmsFullnames
    external :: readclimppt
    external :: read_nldas_climppt

! __________________________________________________________________


   climppt_select = .false.
   do n = 1, LDT_rc%nnest
      if( LDT_climate_struc(n)%climppt%selectOpt == 1 ) then
        climppt_select = .true.
      endif
   enddo
   ! Climate Tmin and Tmax sources not supported at this time.
   ! - Initializing sources to 0 for now.
   climtmin_select = .false.
   LDT_climate_struc%climtmin%selectOpt = 0
   climtmax_select = .false.
   LDT_climate_struc%climtmax%selectOpt = 0

   if( climppt_select ) then
     write(LDT_logunit,*)" - - - - - - - - Climate Downscaling Parameters - - - - - - - - -"
   endif

 !- Read in config file entries:
    if( climppt_select ) then
      call ESMF_ConfigFindLabel(LDT_config,"PPT climatology maps:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,LDT_climate_struc(n)%climpptdir,rc=rc)
         call LDT_verify(rc,'PPT climatology maps: not specified')
      enddo

      do n=1,LDT_rc%nnest
         if (LDT_climate_struc(n)%output_climppt_ratio >0) then
            ! EMK...Legacy NLDAS logic
            if (trim(LDT_climate_struc(n)%climppt2%source) .eq. "NLDAS") then
               call ESMF_ConfigFindLabel(LDT_config,"NLDAS PPT climatology maps:",rc=rc)
               call ESMF_ConfigGetAttribute(LDT_config,LDT_climate_struc(n)%climpptdir2,rc=rc)
               call ESMF_ConfigFindLabel(LDT_config,"NLDAS PPT climatology maps:",rc=rc)
            else
               call ESMF_ConfigFindLabel(LDT_config,"PPT climatology maps 2:",rc=rc)
               call ESMF_ConfigGetAttribute(LDT_config,LDT_climate_struc(n)%climpptdir2,rc=rc)
               call ESMF_ConfigFindLabel(LDT_config,"PPT climatology maps 2:",rc=rc)
            end if
            call ESMF_ConfigFindLabel(LDT_config,"Climate params spatial transform 2:",rc=rc)
            call ESMF_ConfigGetAttribute(LDT_config,LDT_climate_struc(n)%clim_gridtransform2,&
                 rc=rc)
            call LDT_verify(rc,&
                 'Climate params spatial transform 2: option not specified in the config file')
         end if
      enddo

      LDT_rc%monthlyData = .false.
      call ESMF_ConfigFindLabel(LDT_config,"PPT climatology interval:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,LDT_climate_struc(n)%climpptInterval,rc=rc)
         call LDT_verify(rc,'PPT climatology interval: not specified')
         if( LDT_climate_struc(n)%climpptInterval == "monthly" ) then
            LDT_rc%monthlyData(n) = .true.
            LDT_climate_struc(n)%climppt%num_times = 12
         else
            write(LDT_logunit,*) "[ERR] NO OTHER TIME INTERVAL, OTHER THAN 'monthly'"
            write(LDT_logunit,*) "   CURRENTLY EXISTS AS THIS TIME FOR CLIM PPT."
            write(LDT_logunit,*) " Program stopping ..."
            call LDT_endrun
         endif
      enddo
    ! Set units and full names:
      do n=1,LDT_rc%nnest
         LDT_climate_struc(n)%climppt%units="mm"
         call setClimateParmsFullnames( n, "climppt", &
                 LDT_climate_struc(n)%climppt%source )
      enddo
    endif

    if( LDT_climate_struc(1)%climtmin%selectOpt == 1 ) then
      call ESMF_ConfigFindLabel(LDT_config,"TMIN climatology maps:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,LDT_climate_struc(n)%climtmindir,rc=rc)
         call LDT_verify(rc,'TMIN climatology maps: not specified')
      enddo
    endif
    if( LDT_climate_struc(1)%climtmax%selectOpt == 1 ) then
      call ESMF_ConfigFindLabel(LDT_config,"TMAX climatology maps:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,LDT_climate_struc(n)%climtmaxdir,rc=rc)
         call LDT_verify(rc,'TMAX climatology maps: not specified')
      enddo
    end if

    if( climppt_select ) then
      call ESMF_ConfigFindLabel(LDT_config,"Climate params spatial transform:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,LDT_climate_struc(n)%clim_gridtransform,&
              rc=rc)
         call LDT_verify(rc,&
              'Climate params spatial transform: option not specified in the config file')
      enddo
    end if


!- Allocate fields for both LIS running and forcing domains:
   do n = 1, LDT_rc%nnest

   !- Precipitation (PPT) Climatology:
      if( LDT_climate_struc(n)%climppt%selectOpt > 0 ) then 
       !- LIS run domain:
          LDT_climate_struc(n)%climppt%vlevels = LDT_climate_struc(n)%climppt%num_times

          allocate(LDT_climate_struc(n)%climppt%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_climate_struc(n)%climppt%vlevels))

       !- Forcing domain:
          if( LDT_rc%nmetforc > 0 ) then
            do m = 1, LDT_rc%nmetforc_parms
               LDT_force_struc(n,m)%climppt%num_times = LDT_climate_struc(n)%climppt%num_times
               LDT_force_struc(n,m)%climppt%vlevels   = LDT_climate_struc(n)%climppt%num_times
          
               allocate(LDT_force_struc(n,m)%climppt%value(&
                   LDT_rc%met_nc(m), LDT_rc%met_nr(m), &
                   LDT_force_struc(n,m)%climppt%vlevels))       
            end do
          end if
          if ( LDT_climate_struc(n)%output_climppt_ratio>0)then
          allocate(LDT_climate_struc(n)%climppt2%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_climate_struc(n)%climppt%vlevels))
          end if
       endif

#if 0
    !- Max. Temperature Climatology:
       if( LDT_climate_struc(n)%climtmax%selectOpt > 0 ) then 
       !- LIS run domain:
          LDT_climate_struc(n)%climtmax%vlevels = LDT_climate_struc(n)%climtmax%num_times
          LDT_climate_struc(n)%climtmax%num_times = LDT_climate_struc(n)%climtmax%num_times
          allocate(LDT_climate_struc(n)%climtmax%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_climate_struc(n)%climtmax%vlevels))

       !- Forcing domain:
          if( LDT_rc%nmetforc > 0 ) then
            do m = 1, LDT_rc%nmetforc_parms
               LDT_force_struc(n,m)%climtmax%num_times = LDT_climate_struc(n)%climtmax%num_times
               LDT_force_struc(n,m)%climtmax%vlevels = LDT_climate_struc(n)%climtmax%num_times
               allocate(LDT_force_struc(n,m)%climtmax%value(&
                   LDT_rc%met_nc(m), LDT_rc%met_nr(m),&
                   LDT_force_struc(n,m)%climtmax%vlevels))       
            end do
          end if
       endif
#endif

#if 0
    !- Min. Temperature Climatology:
       if( LDT_climate_struc(n)%climtmin%selectOpt > 0 ) then 
       !- LIS run domain:
          LDT_climate_struc(n)%climtmin%vlevels = LDT_climate_struc(n)%climtmin%num_times
          LDT_climate_struc(n)%climtmin%num_times = LDT_climate_struc(n)%climtmin%num_times
          allocate(LDT_climate_struc(n)%climtmin%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_climate_struc(n)%climtmin%vlevels))

       !- Forcing domain:
          if( LDT_rc%nmetforc > 0 ) then
            do m = 1, LDT_rc%nmetforc_parms
               LDT_force_struc(n,m)%climtmin%num_times = LDT_climate_struc(n)%climtmin%num_times
               LDT_force_struc(n,m)%climtmin%vlevels   = LDT_climate_struc(n)%climtmin%num_times
               allocate(LDT_force_struc(n,m)%climtmin%value(&
                   LDT_rc%met_nc(m), LDT_rc%met_nr(m), &
                   LDT_force_struc(n,m)%climtmin%vlevels))       
            end do
          end if
       endif  ! End tmin check
#endif

    enddo  ! end nest loop

 !- Read in data fields:
    do n=1,LDT_rc%nnest

    !- PPT Climatology files:
       if( climppt_select ) then 

!          call LDT_gridOptChecks( n, "Climate PPT", &
!               LDT_climate_struc(n)%clim_gridtransform, &
!               LDT_climate_struc(n)%clim_proj, &
!               clim_gridDesc(n,9) )

          call LDT_climateOptChecks( "Climate PPT", &
               LDT_climate_struc(n)%clim_gridtransform )
          
       !- Monthly (interval) files:
          if( LDT_climate_struc(n)%climpptInterval == "monthly" ) then  
             do k = 1, LDT_climate_struc(n)%climppt%vlevels  ! months
                select case( trim(LDT_climate_struc(n)%climppt%source) )
                  case( "PRISM" ) 
                    LDT_climate_struc(n)%climpptfile = &
                        trim(LDT_climate_struc(n)%climpptdir)//'.'//&
                        trim(months(k))//'.txt'
                  case( "WORLDCLIM" ) 
                    LDT_climate_struc(n)%climpptfile = &
                        trim(LDT_climate_struc(n)%climpptdir)//&
                        trim(mon2d(k))
                 case ("NAFPA_BACK_GALWEM")
                    LDT_climate_struc(n)%climpptfile = &
                         trim(LDT_climate_struc(n)%climpptdir)
                    read(mon2d(k),'(I3)') LDT_climate_struc(n)%climpptimonth
!                 case ("NAFPA_BACK_GFS")
!                    LDT_climate_struc(n)%climpptfile = &
!                         trim(LDT_climate_struc(n)%climpptdir)
!                    read(mon2d(k),'(I3)') LDT_climate_struc(n)%climpptimonth
                 case default
                    write(LDT_logunit,*) "[ERR] PPT Climatology Source Not Recognized"
                    write(LDT_logunit,*) trim(LDT_climate_struc(n)%climppt%source)
                    write(LDT_logunit,*) " Program stopping ..."
                    call LDT_endrun
                end select

                if (LDT_climate_struc(n)%output_climppt_ratio>0) then

                   select case( trim(LDT_climate_struc(n)%climppt2%source) )
                   case( "NLDAS" )
                      LDT_climate_struc(n)%climpptfile2 = &
                           trim(LDT_climate_struc(n)%climpptdir2)//'.'//&
                           trim(months(k))//'.txt'
                   !case ("NAFPA_BACK_GALWEM")
                   !   LDT_climate_struc(n)%climpptfile2 = &
                   !        trim(LDT_climate_struc(n)%climpptdir2)
                   !   read(mon2d(k),'(I3)') LDT_climate_struc(n)%climpptimonth
                   case ("NAFPA_BACK_GFS")
                      LDT_climate_struc(n)%climpptfile2 = &
                           trim(LDT_climate_struc(n)%climpptdir2)
                      read(mon2d(k),'(I3)') LDT_climate_struc(n)%climpptimonth
                   case default
                      write(LDT_logunit,*) "[ERR] PPT Climatology Source 2 Not Recognized"
                      write(LDT_logunit,*) trim(LDT_climate_struc(n)%climppt2%source)
                      write(LDT_logunit,*) " Program stopping ..."
                      call LDT_endrun
                   end select


                end if

               !- Read monthly climatology files:
               write(LDT_logunit,*) "Reading "//trim(LDT_climate_struc(n)%climpptfile)

             !- For LIS Run domain:
               call readclimppt( trim(LDT_climate_struc(n)%climppt%source)//char(0),&
                     n, LDT_rc%lnc(n), LDT_rc%lnr(n), LDT_rc%gridDesc(n,:), &
                     LDT_climate_struc(n)%climppt%value(:,:,k) )

                !if the second climatology file exists, calculuate the ratios of the two
                if (LDT_climate_struc(n)%output_climppt_ratio>0) then

                   !EMK Preserve legacy NLDAS logic
                   if (trim(LDT_climate_struc(n)%climppt2%source) .eq. "NLDAS") then
                      !- Read NLDAS PPT climatology for LIS Run domain:
                      write(LDT_logunit,*) "Reading "//trim(LDT_climate_struc(n)%climpptfile2)
                      call read_NLDAS_climppt(n, LDT_rc%lnc(n), LDT_rc%lnr(n), LDT_rc%gridDesc(n,:), &
                           LDT_climate_struc(n)%climppt2%value(:,:,k))
                   else
                      call readclimppt( trim(LDT_climate_struc(n)%climppt2%source)//char(0),&
                           n, LDT_rc%lnc(n), LDT_rc%lnr(n), LDT_rc%gridDesc(n,:), &
                           LDT_climate_struc(n)%climppt2%value(:,:,k) )
                   end if

                   !calculate the ratio
                   do ir =1,LDT_rc%lnr(n)
                      do ic=1,LDT_rc%lnc(n)
                         ! Handle data voids
                         if (LDT_climate_struc(n)%climppt%value(ic,ir,k).eq.LDT_rc%udef .or. &
                              LDT_climate_struc(n)%climppt2%value(ic,ir,k).eq.LDT_rc%udef) then
                            ratio = 1.0

                         ! If both datasets average less than a monthly trace,
                         ! treat as a data void.
                         else if (LDT_climate_struc(n)%climppt%value(ic,ir,k) .lt. 0.1 .and. &
                              LDT_climate_struc(n)%climppt%value(ic,ir,k) .lt. 0.1) then
                            ratio = 1.0
                         ! If top dataset has less than a monthly trace, but
                         ! the bottom dataset has more, assume a trace for
                         ! the top (prevent zero bias ratio).
                         else if (LDT_climate_struc(n)%climppt%value(ic,ir,k) .lt. 0.1) then
                            ratio = 0.1 / &
                                 LDT_climate_struc(n)%climppt2%value(ic,ir,k) + 0.1
                         ! If bottom dataset has less than a monthly trace, but
                         ! the top dataset has more, assume a trace for
                         ! the bottom (prevent divide by zero).
                         else if (LDT_climate_struc(n)%climppt2%value(ic,ir,k) .lt. 0.1) then
                            ratio = LDT_climate_struc(n)%climppt%value(ic,ir,k) / &
                                 (0.1)
                         ! Normal bias ratio.
                         else
                            ratio = LDT_climate_struc(n)%climppt%value(ic,ir,k) / &
                                 LDT_climate_struc(n)%climppt2%value(ic,ir,k)

                         end if

                         ! if (ratio .gt. 10) then
                         !    write(LDT_logunit,*) &
                         !         'EMK: ic,ir,k,top,bottom,diff,ratio = ', ic, ir, k, &
                         !         LDT_climate_struc(n)%climppt%value(ic,ir,k), &
                         !         LDT_climate_struc(n)%climppt2%value(ic,ir,k), &
                         !         LDT_climate_struc(n)%climppt%value(ic,ir,k) - &
                         !         LDT_climate_struc(n)%climppt2%value(ic,ir,k), &
                         !         ratio
                         ! end if
                         LDT_climate_struc(n)%climppt%value(ic,ir,k) = ratio
                      end do
                   end do
                end if

             !- For Forcing domains:
                if( LDT_rc%nmetforc > 0 ) then
                  do m = 1, LDT_rc%nmetforc_parms

                  !- Transfer common parameter attribute entries:
                  !  from LIS-run domain to MET-run domains.
                  ! ( input_parmattribs -> output_parmattribs ) 
                     call populate_param_attribs( &
!                          "PPTCLIM_"//trim(LDT_rc%metforc_parms(m)),   &
                          "PPTCLIM_"//trim(LDT_rc%metforc_parmsrc(m)),   &
                           LDT_climate_struc(n)%climppt%standard_name//&
!                               " for "//trim(LDT_rc%metforc_parms(m)), "mm", &
                               " for "//trim(LDT_rc%metforc_parmsrc(m)), "mm", &
                           LDT_climate_struc(n)%climppt,               &
                           LDT_force_struc(n,m)%climppt )


                     write(LDT_logunit,*) "Processing "//trim(LDT_climate_struc(n)%climpptfile)
                     write(LDT_logunit,*) " for metforcing grid: "//trim(LDT_rc%metforc_parms(m))

                     call readclimppt( trim(LDT_climate_struc(n)%climppt%source)//char(0),&
                          n, LDT_rc%met_nc(m), LDT_rc%met_nr(m), LDT_rc%met_gridDesc(m,:), &
                          LDT_force_struc(n,m)%climppt%value(:,:,k) )
                  end do
                end if
                write(LDT_logunit,*) "Done reading "//trim(LDT_climate_struc(n)%climpptfile)
             enddo
             
          endif  ! End monthly check
       endif     ! End ClimPPT check

    !- Min Temperature Climatology Files:
!       if( LDT_climate_struc(n)%climtmin%selectOpt == 1 ) then 
!       endif
    !- Max Temperature Climatology Files:
!       if( LDT_climate_struc(n)%climtmax%selectOpt == 1 ) then 
!       endif

    enddo

  end subroutine LDT_climateparms_init

  subroutine LDT_climateparms_writeHeader(n,ftn,dimID,met_dimID,monthID)

    integer     :: n 
    integer     :: ftn
    integer     :: dimID(3)
    integer     :: met_dimID(LDT_rc%nmetforc_parms,3)
    integer     :: monthID

    integer     :: m
    integer     :: t_dimID(3)

! - Climatological PPT Fields:
    if( LDT_climate_struc(n)%climppt%selectOpt.gt.0 ) then 

       if( LDT_climate_struc(n)%climpptInterval.eq."monthly" ) then  ! monthly
          t_dimID(3) = monthID       
       endif
    !- LIS run domain:
       t_dimID(1) = dimID(1)
       t_dimID(2) = dimID(2)
       call LDT_writeNETCDFdataHeader(n,ftn,t_dimID,&
            LDT_climate_struc(n)%climppt)

    !- Forcing domain:
       if( LDT_rc%nmetforc > 0 .and. LDT_climate_struc(n)%output_climppt_ratio==0) then
         do m = 1, LDT_rc%nmetforc_parms
            t_dimID(1) = met_dimID(m,1)
            t_dimID(2) = met_dimID(m,2)
            call LDT_writeNETCDFdataHeader(n,ftn,t_dimID,&
                 LDT_force_struc(n,m)%climppt)
         end do
       end if
    endif

! - Climatological TMIN Fields:
!          call LDT_writeNETCDFdataHeader(n,ftn,t_dimID,&
!               LDT_climate_struc(n)%climtmin)


!          call LDT_writeNETCDFdataHeader(n,ftn,t_dimID,&
!               LDT_climate_struc(n)%climtmax)


  end subroutine LDT_climateparms_writeHeader

  subroutine LDT_climateparms_writeData(n,ftn)

    integer  :: n 
    integer  :: ftn
    integer  :: m

    if( LDT_climate_struc(n)%climppt%selectOpt.gt.0 ) then 
       call LDT_writeNETCDFdata(n, ftn, LDT_climate_struc(n)%climppt)

    !- Forcing domain:
       if( LDT_rc%nmetforc > 0 .and. LDT_climate_struc(n)%output_climppt_ratio==0 ) then
          do m = 1, LDT_rc%nmetforc_parms
            call LDT_writeNETCDFdata( n, ftn, LDT_force_struc(n,m)%climppt, &
                 LDT_rc%met_nc(m), LDT_rc%met_nr(m) )
         end do
       end if
    endif

!       call LDT_writeNETCDFdata(n,ftn,LDT_force_struc(n,m)%climtmin)
!       call LDT_writeNETCDFdata(n,ftn,LDT_force_struc(n,m)%climtmax)

  end subroutine LDT_climateparms_writeData

end module LDT_climateParmsMod
