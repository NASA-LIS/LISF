!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
module LDT_laisaiMod
!BOP
!
! !MODULE: LDT_laisaiMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read leaf area index (LAI)
!  and stem area index (SAI) data. 
!
!  \subsubsection{Overview}
!  This routines in this module provides routines to read the 
!  LAI/SAI climatology data and allows the users to 
!  specify the frequency of climatology (in months). 
!  The climatological data is temporally interpolated  
!  between months to the current simulation date. 
!
! !REVISION HISTORY:
!
!  08 Aug 2005: Sujay Kumar; Initial implementation
!  12 Feb 2013: KR Arsenault: Modified for CLSM LAI data.
!
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use ESMF
  use LDT_coreMod
  use LDT_historyMod
  use LDT_paramDataMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_paramMaskCheckMod

  implicit none

  PRIVATE

!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LDT_laisai_readParamSpecs
  public :: LDT_laisai_init         ! allocates memory for required structures
  public :: LDT_laisai_writeHeader
  public :: LDT_laisai_writeData
!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LDT_laisai_struc

  type, public :: laisai_type_dec
     real                   :: laisai_gridDesc(20)
     character*50           :: laisai_proj
     character*50           :: laisai_gridtransform
     character(len=LDT_CONST_PATH_LEN)          :: laidir
     character(len=LDT_CONST_PATH_LEN)          :: saidir
     character(len=LDT_CONST_PATH_LEN)          :: laifile
     character(len=LDT_CONST_PATH_LEN)          :: saifile
     character*20           :: laisaiInterval

     character(len=LDT_CONST_PATH_LEN)          :: laimaxfile
     character(len=LDT_CONST_PATH_LEN)          :: laiminfile

     type(LDT_paramEntry) :: lai         ! Leaf area-index (LAI)
     type(LDT_paramEntry) :: laimin      ! Min. LAI (@pixel)
     type(LDT_paramEntry) :: laimax      ! Max. LAI (@pixel)
     type(LDT_paramEntry) :: sai         ! Stem area-index (SAI)

  end type laisai_type_dec

  type(laisai_type_dec), allocatable :: LDT_laisai_struc(:)

contains

  subroutine LDT_laisai_readParamSpecs
    
    character*100     :: source
    integer           :: rc
    integer           :: n
    
    allocate(LDT_laisai_struc(LDT_rc%nnest))
    
    call ESMF_ConfigFindLabel(LDT_config,"LAI data source:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
!       call LDT_warning(rc,"LAI data source: not defined")
       call LDT_set_param_attribs(rc,LDT_laisai_struc(n)%lai,&
            "LAI",source)
       call LDT_set_param_attribs(rc,LDT_laisai_struc(n)%laimin,&
            "LAIMIN",source)
       call LDT_set_param_attribs(rc,LDT_laisai_struc(n)%laimax,&
            "LAIMAX",source)
    enddo

  ! LSM-required parameter check:
    if( index(LDT_rc%lsm,"CLM")  == 1 .or. &
        index(LDT_rc%lsm,"CLSM") == 1 ) then
      if( rc /= 0 ) then
         call LDT_warning(rc,"WARNING: LAI data source: not defined")
      endif
    endif

    call ESMF_ConfigFindLabel(LDT_config,"SAI data source:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
       call LDT_set_param_attribs(rc,LDT_laisai_struc(n)%sai,&
            "SAI",source)
    enddo

  ! LSM-required parameter check:
    if( index(LDT_rc%lsm,"CLM")  == 1 .or. &
        index(LDT_rc%lsm,"CLSM") == 1 ) then
      if( rc /= 0 ) then
         call LDT_warning(rc,"WARNING: SAI data source: not defined")
      endif
    endif

  end subroutine LDT_laisai_readParamSpecs
!BOP
! 
! !ROUTINE: LDT_laisai_init
! \label{LDT_laisai_init}
! 
! !INTERFACE:
  subroutine LDT_laisai_init

! !USES:
    use LDT_fileIOMod, only : LDT_readDomainConfigSpecs
    use LDT_logMod,    only : LDT_verify
    use LDT_paramOptCheckMod, only: LDT_laisaiOptChecks, &
                       LDT_gridOptChecks
 
! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! the LAI/SAI datasets.
! 
!  The routines invoked are: 
!  \begin{description}
!  \end{description}
!
!EOP
   implicit none
   integer   :: n
   integer   :: k, c, r
   integer   :: rc
   character*3 :: months(12)
   data months /'jan','feb','mar','apr','may','jun','jul','aug',&
                'sep','oct','nov','dec'/

   type(LDT_fillopts) :: lai
   type(LDT_fillopts) :: laimax
   type(LDT_fillopts) :: laimin
   type(LDT_fillopts) :: sai
   logical            :: lai_select, check_data1, check_data2
   logical            :: calc_minmaxlai

   real, allocatable           :: laisai_gridDesc(:,:)
   character*50                :: laisai_proj
   character*50,  allocatable  :: laisai_gridtransform(:)
   character*100, allocatable  :: laidir(:)
   character*140, allocatable  :: laifile(:)
   character*100, allocatable  :: saidir(:)
   character*140, allocatable  :: saifile(:)
   character*20,  allocatable  :: laisaiInterval(:)

! ___________________________________________________________________

    lai_select = .false.
    do n=1,LDT_rc%nnest
       if(LDT_laisai_struc(n)%lai%selectOpt.gt.0) then
          lai_select = .true.
       endif
    enddo
    if(lai_select) then
       write(LDT_logunit,*)" - - - - - - - - - - LAI/SAI Parameters - - - - - - - - - -"
       call setlaiattribs(trim(LDT_laisai_struc(1)%lai%source)//char(0))
    endif

   allocate(laisai_gridDesc(LDT_rc%nnest,20))
   allocate(laisai_gridtransform(LDT_rc%nnest))         
   allocate(laifile(LDT_rc%nnest))
   allocate(laidir(LDT_rc%nnest))
   allocate(saifile(LDT_rc%nnest))
   allocate(saidir(LDT_rc%nnest))
   allocate(laisaiInterval(LDT_rc%nnest))

    do n=1,LDT_rc%nnest

       LDT_laisai_struc(n)%lai%vlevels = LDT_laisai_struc(n)%lai%num_times
       LDT_laisai_struc(n)%sai%vlevels = LDT_laisai_struc(n)%lai%num_times ! Set to lai for now

       LDT_laisai_struc(n)%laimax%vlevels = LDT_laisai_struc(n)%laimax%num_times
       LDT_laisai_struc(n)%laimin%vlevels = LDT_laisai_struc(n)%laimin%num_times

       if( LDT_laisai_struc(n)%lai%selectOpt.eq.1 ) then 
          allocate(LDT_laisai_struc(n)%lai%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_laisai_struc(n)%lai%vlevels))       
       endif
       if( LDT_laisai_struc(n)%sai%selectOpt.eq.1 ) then 
          allocate(LDT_laisai_struc(n)%sai%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_laisai_struc(n)%sai%vlevels))       
       endif
       if(LDT_laisai_struc(n)%laimax%selectOpt.gt.0) then 
          allocate(LDT_laisai_struc(n)%laimax%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_laisai_struc(n)%laimax%vlevels))       
       endif
       if(LDT_laisai_struc(n)%laimin%selectOpt.gt.0) then 
          allocate(LDT_laisai_struc(n)%laimin%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_laisai_struc(n)%laimin%vlevels))       
       endif
    end do 

 !- LAI:
    if( lai_select ) then 
       call ESMF_ConfigFindLabel(LDT_config,"LAI map:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,laidir(n),rc=rc)
          call LDT_verify(rc,'LAI map: not specified')
          LDT_laisai_struc(n)%laidir = laidir(n)
       enddo
     ! Set units and full names:
       do n=1,LDT_rc%nnest
          LDT_laisai_struc(n)%lai%units="-"
          call setLAISAIParmsFullnames( n, "lai", &
                  LDT_laisai_struc(n)%lai%source )
       enddo
    end if
 !- SAI:
    if( LDT_laisai_struc(1)%SAI%selectOpt.eq.1 ) then 
       call ESMF_ConfigFindLabel(LDT_config,"SAI map:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,saidir(n),rc=rc)
          call LDT_verify(rc,'SAI map: not specified')
          LDT_laisai_struc(n)%saidir = saidir(n)
       enddo
     ! Set units and full names:
       do n=1,LDT_rc%nnest
          LDT_laisai_struc(n)%sai%units="-"
          call setLAISAIParmsFullnames( n, "sai", &
                  LDT_laisai_struc(n)%sai%source )
       enddo
    endif

 !- Determine additional LAI/SAI config file entries:
    if( lai_select ) then
       call ESMF_ConfigFindLabel(LDT_config,"LAI/SAI climatology interval:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,laisaiInterval(n),rc=rc)
          call LDT_verify(rc,'LAI/SAI climatology interval: not specified')

         if( trim(laisaiInterval(n)) .ne. "monthly" ) then
           write(LDT_logunit,*) "ERR: 'monthly' LAI interval option specified."
           write(LDT_logunit,*) "     Please select 'monthly' at this time. ... Stopping."
           call LDT_endrun
         endif
         if( laisaiInterval(n) == "monthly" .and. &
             LDT_laisai_struc(n)%lai%vlevels .ne. 12 )then
            write(LDT_logunit,*) "ERR: The 'monthly' LAI interval option should have '12' in the"
            write(LDT_logunit,*) "     parameter attribs table.  Please change to '12' there.  Stopping."
            call LDT_endrun
         endif
         LDT_laisai_struc(n)%laisaiInterval = laisaiInterval(n)
       enddo

       lai%filltype = "none"
       call ESMF_ConfigGetAttribute(LDT_config, lai%filltype, &
            label="LAI/SAI fill option:",rc=rc)
       call LDT_verify(rc,"LAI/SAI fill option: option not specified in the config file")

       if( lai%filltype == "average" .or. lai%filltype == "neighbor" ) then
          call ESMF_ConfigGetAttribute(LDT_config, lai%fillradius, &
               label="LAI/SAI fill radius:",rc=rc)
          call LDT_verify(rc,"LAI/SAI fill radius: option not specified in the config file")

          if( lai_select ) then
             call ESMF_ConfigGetAttribute(LDT_config, lai%fillvalue, &
                  label="LAI fill value:",rc=rc)
             call LDT_verify(rc,"LAI fill value: option not specified in the config file")
          endif
          if( LDT_laisai_struc(1)%laimax%selectOpt == 1 ) then
             call ESMF_ConfigGetAttribute(LDT_config, laimax%fillvalue, &
                  label="LAI maximum fill value:",rc=rc)
             call LDT_verify(rc,"LAI maximum fill value: option not specified in the config file")
          endif
          if( LDT_laisai_struc(1)%laimin%selectOpt == 1 ) then
             call ESMF_ConfigGetAttribute(LDT_config, laimin%fillvalue, &
                  label="LAI minimum fill value:",rc=rc)
             call LDT_verify(rc,"LAI minimum fill value: option not specified in the config file")
          endif
          if( LDT_laisai_struc(1)%sai%selectOpt == 1 ) then
            call ESMF_ConfigGetAttribute(LDT_config, sai%fillvalue, &
                 label="SAI fill value:",rc=rc)
            call LDT_verify(rc,"SAI fill value: option not specified in the config file")
          endif
       elseif( lai%filltype == "none" ) then
         write(LDT_logunit,*) "[INFO] 'NONE' Parameter-Mask Agreement Option Selected for LAI/SAI"
       else
         write(LDT_logunit,*) "[ERR] Fill option for LAI/SAI is not valid: ",trim(lai%filltype)
         write(LDT_logunit,*) "  Please select one of these:  none, neighbor or average "
         write(LDT_logunit,*) "  Programming stopping ..."
         call LDT_endrun
       end if

       call ESMF_ConfigGetAttribute(LDT_config,laisai_proj,&
            label="LAI/SAI map projection:",rc=rc)
       call LDT_verify(rc,'LAI/SAI projection: option not specified in the config file')
       LDT_laisai_struc(:)%laisai_proj = laisai_proj

       call ESMF_ConfigFindLabel(LDT_config,"LAI/SAI spatial transform:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,laisai_gridtransform(n),&
               rc=rc)
          call LDT_verify(rc,'LAI/SAI spatial transform: option not specified in the config file')
          LDT_laisai_struc(n)%laisai_gridtransform = laisai_gridtransform(n)
       enddo

       call LDT_readDomainConfigSpecs("LAI/SAI",laisai_proj,&
            laisai_gridDesc)
    endif
    do n=1,LDT_rc%nnest
       LDT_laisai_struc(n)%laisai_gridDesc = laisai_gridDesc(n,:)
    enddo

!- Check if LAI min and max fields are selected:
    check_data1 = .false.
    check_data2 = .false.
    do n=1,LDT_rc%nnest
       if(LDT_laisai_struc(n)%laimax%selectOpt.gt.0 ) then
          check_data1 = .true.
       endif
       if(LDT_laisai_struc(n)%laimin%selectOpt.gt.0 ) then
          check_data2 = .true.
       endif
    enddo
    
!- Decide to either calculate or read-in min/max LAI values:
    calc_minmaxlai = .false.
    if( check_data1 .or. check_data2 ) then
       call ESMF_ConfigGetAttribute(LDT_config,calc_minmaxlai,&
            label="Calculate min-max LAI:",rc=rc)
       call LDT_verify(rc,'Calculate min-max LAI: not specified')
       if( .NOT. calc_minmaxlai ) then
          if( check_data1 ) then
             !- Max LAI file:
             call ESMF_ConfigFindLabel(LDT_config,"LAI maximum map:",rc=rc)
             do n=1,LDT_rc%nnest
                call ESMF_ConfigGetAttribute(LDT_config,&
                     LDT_laisai_struc(n)%laimaxfile,rc=rc)
                call LDT_verify(rc,'LAI maximum map: not specified')
             enddo
          endif
          if( check_data2 ) then
             !- Min LAI file:
             call ESMF_ConfigFindLabel(LDT_config,"LAI minimum map:",rc=rc)
             do n=1,LDT_rc%nnest
                call ESMF_ConfigGetAttribute(LDT_config,&
                LDT_laisai_struc(n)%laiminfile,rc=rc)
                call LDT_verify(rc,'LAI minimum map: not specified')
             enddo
          endif
       endif
     ! Set units and full names:
       do n=1,LDT_rc%nnest
          LDT_laisai_struc(n)%laimin%units="-"
          LDT_laisai_struc(n)%laimax%units="-"
          call setLAISAIParmsFullnames( n, "laimin", &
                  LDT_laisai_struc(n)%lai%source )
          call setLAISAIParmsFullnames( n, "laimax", &
                  LDT_laisai_struc(n)%lai%source )
       enddo
    endif
    
    do n=1,LDT_rc%nnest

    !- Check LAI grid options:
       if(( LDT_laisai_struc(n)%lai%selectOpt == 1 )  .or. & 
          ( LDT_laisai_struc(n)%sai%selectOpt == 1 )) then 
          call LDT_gridOptChecks( n, "LAI/SAI", &
               laisai_gridtransform(n), laisai_proj, &
               laisai_gridDesc(n,9) )
          
          call LDT_laisaiOptChecks( "LAI/SAI", laisai_proj,&
               laisai_gridtransform(n) )
       endif

    !- Read in LAI/SAI files:
       if(LDT_laisai_struc(n)%lai%selectOpt.eq.1) then 

          if(laisaiInterval(n).eq."monthly") then !monthly
             LDT_rc%monthlyData(n) = .true.

          !- Read single-file monthly clim LAI: 
             if( trim(LDT_laisai_struc(n)%lai%source) == "CLSMF2.5" ) then
                LDT_laisai_struc(n)%laifile = trim(laidir(n))
                write(LDT_logunit,*) "Reading single-file, monthly climatologies for: "&
                     //trim(LDT_laisai_struc(n)%laifile)
                !call readlai( trim(LDT_laisai_struc(n)%lai%source)//char(0),&
                !     n, LDT_laisai_struc(n)%lai%value, &
                !     LDT_LSMparam_struc(n)%landmask%value )
                call readlai( trim(LDT_laisai_struc(n)%lai%source)//char(0),&
                     n, LDT_laisai_struc(n)%lai%value)
                write(LDT_logunit,*) "Done reading file - "//&
                     trim(LDT_laisai_struc(n)%laifile)
                
          !- Read multi-file monthly LAI: 
             else  ! Other LAI sources
                do k=1,LDT_laisai_struc(n)%lai%vlevels
                   LDT_laisai_struc(n)%laifile = trim(laidir(n))//'.'//&
                       trim(months(k))//'.1gd4r'
                   write(LDT_logunit,*) 'Reading '//trim(LDT_laisai_struc(n)%laifile)
                   call readlai(trim(LDT_laisai_struc(n)%lai%source)//char(0),&
                        n,LDT_laisai_struc(n)%lai%value(:,:,k))
                   write(LDT_logunit,*) 'Done reading '//&
                        trim(LDT_laisai_struc(n)%laifile)
                enddo
             end if
          endif  ! End interval check

          if( lai%filltype == "average" .or. lai%filltype == "neighbor" ) then
             write(LDT_logunit,*) "Checking/filling mask values for: ", &
                                trim(LDT_laisai_struc(n)%lai%short_name)
             write(fill_logunit,*) "Checking/filling mask values for: ", &
                                trim(LDT_laisai_struc(n)%lai%short_name)
             lai%watervalue = LDT_rc%udef
!       fill_option = "average"
!       fill_value = 1.0
!       fill_rad = 2.
             call LDT_contIndivParam_Fill(n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
                 laisai_gridtransform(n), LDT_laisai_struc(n)%lai%num_times, &
                 LDT_laisai_struc(n)%lai%value, lai%watervalue,       &
                 LDT_LSMparam_struc(n)%landmask2%value, lai%filltype,      &
                 lai%fillvalue, lai%fillradius )
          endif
       endif   ! end LAI read
       
    !- Maximum LAI file:
       if( LDT_laisai_struc(n)%laimax%selectOpt == 1 ) then 

       !- Calculate maximum LAI values:
          if( calc_minmaxlai ) then
            do r = 1, LDT_rc%lnr(n)
               do c = 1, LDT_rc%lnc(n)
                  LDT_laisai_struc(n)%laimax%value(c,r,1) = &
                      maxval(LDT_laisai_struc(n)%lai%value(c,r,:))
               enddo
            enddo
       !- Read maximum LAI file:
          else
             call readlaimax(trim(LDT_laisai_struc(n)%lai%source)//char(0),&
                             n,LDT_laisai_struc(n)%laimax%value(:,:,1))
          endif 
       !- Perform mask-parameter consistency "fill" options:
          if( lai%filltype == "average" .or. lai%filltype == "neighbor" ) then
             write(LDT_logunit,*) "Checking/filling mask values for: ", &
                  trim(LDT_laisai_struc(n)%laimax%short_name)
             write(fill_logunit,*) "Checking/filling mask values for: ", &
                  trim(LDT_laisai_struc(n)%laimax%short_name)
             lai%watervalue = LDT_rc%udef
             call LDT_contIndivParam_Fill(n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
                  LDT_laisai_struc(n)%laisai_gridtransform, &
                  LDT_laisai_struc(n)%laimax%num_times, &
                  LDT_laisai_struc(n)%laimax%value, lai%watervalue,      &
                  LDT_LSMparam_struc(n)%landmask2%value, lai%filltype,      &
                  laimax%fillvalue, lai%fillradius )
          endif
       endif

    !- Minimum LAI Values:
       if( LDT_laisai_struc(n)%laimin%selectOpt == 1 ) then 
       !- Calculate min LAI values:
          if( calc_minmaxlai ) then
            do r = 1, LDT_rc%lnr(n)
               do c = 1, LDT_rc%lnc(n)
                  LDT_laisai_struc(n)%laimin%value(c,r,1) = &
                      minval(LDT_laisai_struc(n)%lai%value(c,r,:))
               enddo
            enddo
       !- Read minimum LAI file:
          else
            call readlaimin(trim(LDT_laisai_struc(n)%lai%source)//char(0),&
                            n,LDT_laisai_struc(n)%laimin%value(:,:,1))
          endif
          
       !- Perform mask-parameter consistency "fill" options:
          if( lai%filltype == "average" .or. lai%filltype == "neighbor" ) then
             write(LDT_logunit,*) "Checking/filling mask values for: ", &
                  trim(LDT_laisai_struc(n)%laimin%short_name)
             write(fill_logunit,*) "Checking/filling mask values for: ", &
                  trim(LDT_laisai_struc(n)%laimin%short_name)
             lai%watervalue = LDT_rc%udef
             call LDT_contIndivParam_Fill(n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
                  LDT_laisai_struc(n)%laisai_gridtransform, &
                  LDT_laisai_struc(n)%laimin%num_times, &
                  LDT_laisai_struc(n)%laimin%value, lai%watervalue,      &
                  LDT_LSMparam_struc(n)%landmask2%value, lai%filltype,      &
                  laimin%fillvalue, lai%fillradius )
          endif
       endif
       if( LDT_laisai_struc(n)%sai%selectOpt.eq.1 ) then 
          if( laisaiInterval(n).eq."monthly" ) then   ! monthly
             LDT_rc%monthlyData(n) = .true.
             do k = 1, LDT_laisai_struc(n)%sai%vlevels
                LDT_laisai_struc(n)%saifile = trim(saidir(n))//'.'//&
                     trim(months(k))//'.1gd4r'
                write(LDT_logunit,*) 'Reading '//trim(LDT_laisai_struc(n)%saifile)
                call readsai(trim(LDT_laisai_struc(n)%sai%source)//char(0),&
                     n,LDT_laisai_struc(n)%sai%value(:,:,k))
                write(LDT_logunit,*) 'Done reading '//trim(LDT_laisai_struc(n)%saifile)
             enddo
          endif   ! end interval condition

          if( lai%filltype == "average" .or. lai%filltype == "neighbor" ) then
             write(LDT_logunit,*) "Checking/filling mask values for: ", &
                                trim(LDT_laisai_struc(n)%sai%short_name)
             write(fill_logunit,*) "Checking/filling mask values for: ", &
                                trim(LDT_laisai_struc(n)%sai%short_name)
             lai%watervalue = LDT_rc%udef
!       fill_option = "average"
!       fill_value = 1.0
!       fill_rad = 2.
             call LDT_contIndivParam_Fill(n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
                  LDT_laisai_struc(n)%laisai_gridtransform, &
                  LDT_laisai_struc(n)%sai%num_times, &
                  LDT_laisai_struc(n)%sai%value, lai%watervalue,         &
                  LDT_LSMparam_struc(n)%landmask2%value, lai%filltype,      &
                  sai%fillvalue, lai%fillradius )
          endif

       endif ! end SAI read

    enddo    ! end nest loop

  end subroutine LDT_laisai_init
  
  subroutine LDT_laisai_writeHeader(n,ftn,dimID,monthID)

    integer            :: n 
    integer            :: ftn
    integer            :: dimID(3)
    integer            :: monthID

    integer            :: t_dimID(3)

    if(LDT_laisai_struc(n)%lai%selectOpt.eq.1.or.&
       LDT_laisai_struc(n)%sai%selectOpt.eq.1) then 

       if(LDT_laisai_struc(n)%laisaiInterval.eq."monthly") then !monthly
          t_dimID(1) = dimID(1)
          t_dimID(2) = dimID(2)
          t_dimID(3) = monthID       
          call LDT_writeNETCDFdataHeader(n,ftn,t_dimID,&
               LDT_laisai_struc(n)%lai)
          call LDT_writeNETCDFdataHeader(n,ftn,t_dimID,&
               LDT_laisai_struc(n)%laimin)
          call LDT_writeNETCDFdataHeader(n,ftn,t_dimID,&
               LDT_laisai_struc(n)%laimax)
          call LDT_writeNETCDFdataHeader(n,ftn,t_dimID,&
               LDT_laisai_struc(n)%sai)

       endif

       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"LAISAI_DATA_INTERVAL", &
            LDT_laisai_struc(n)%laisaiInterval))
    endif

  end subroutine LDT_laisai_writeHeader

  subroutine LDT_laisai_writeData(n,ftn)


    integer          :: n 
    integer          :: ftn

    if(LDT_laisai_struc(n)%LAI%selectOpt.eq.1.or.&
         LDT_laisai_struc(n)%sai%selectOpt.eq.1) then 
       call LDT_writeNETCDFdata(n,ftn,LDT_laisai_struc(n)%lai)
       call LDT_writeNETCDFdata(n,ftn,LDT_laisai_struc(n)%laimin)
       call LDT_writeNETCDFdata(n,ftn,LDT_laisai_struc(n)%laimax)
       call LDT_writeNETCDFdata(n,ftn,LDT_laisai_struc(n)%sai)
    endif

  end subroutine LDT_laisai_writeData

end module LDT_laisaiMod
