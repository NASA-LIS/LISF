!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module Noah_parmsMod
!BOP
!
! !MODULE: Noah_parmsMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read greenness fraction
!  data. 
!  \subsubsection{Overview}
!  This routines in this module provides routines to read the 
!  greenness fraction climatology data and allows the users to 
!  specify the frequency of climatology (in months). 
!  The climatological data is temporally interpolated  
!  between months to the current simulation date. 
!
! !REVISION HISTORY:
!
!  08 Aug 2005: Sujay Kumar; Initial implementation
!  04 Aug 2012: K. Arsenault: Made updates to Tbot inputs  
!  27 Aug 2021: Sarith Mahanama: MMF groundwater parameters were added.
  
  use ESMF
  use LDT_coreMod
  use LDT_historyMod
  use LDT_paramDataMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_paramMaskCheckMod
  use MMF_groundwater, ONLY : MMF_BCsReader, cell_area

  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: NoahParms_init    !allocates memory for required structures
  public :: NoahParms_writeHeader
  public :: NoahParms_writeData

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: Noah_struc

  type, public :: noah_type_dec

     real           :: tbot_gridDesc(20)
     character(len=LDT_CONST_PATH_LEN)  :: tbotfile
     character*50   :: tbot_gridtransform
     character*50   :: tbot_topocorr

     real           :: slopetype_gridDesc(20)
     character(len=LDT_CONST_PATH_LEN)  :: slopetypefile
     character*50   :: slopetype_gridtransform

     character*50   :: tbot_proj
     character*50   :: slopetype_proj

     character(len=LDT_CONST_PATH_LEN)  :: mmf_fdepth_dir
     character(len=LDT_CONST_PATH_LEN)  :: mmf_rechclim_dir
     character(len=LDT_CONST_PATH_LEN)  :: mmf_riverbed_dir
     character(len=LDT_CONST_PATH_LEN)  :: mmf_eqwtd_dir
     character(len=LDT_CONST_PATH_LEN)  :: mmf_hgtm_dir
     character*50   :: mmf_transform
      
     real           :: pblh_value

     ! -  Noah LSM-specific:
     type(LDT_paramEntry) :: tbot        ! Bottom temperature (Noah)
     type(LDT_paramEntry) :: slopetype   ! Slope type index (Noah)
     type(LDT_paramEntry) :: pblh        ! Planetary Boundary Layer Height (Noah-MP)
     type(LDT_paramEntry) :: fdepth
     type(LDT_paramEntry) :: rechclim
     type(LDT_paramEntry) :: riverbed
     type(LDT_paramEntry) :: eqwtd
     type(LDT_paramEntry) :: areaxy
     type(LDT_paramEntry) :: hgtm
     
  end type noah_type_dec

  type(noah_type_dec), allocatable :: Noah_struc(:)

contains

!BOP
! 
! !ROUTINE: NoahParms_init
! \label{NoahParms_init}
! 
! !INTERFACE:

  subroutine NoahParms_init (flag)
    ! !USES:
    use LDT_fileIOMod, only : LDT_readDomainConfigSpecs
    use LDT_logMod,    only : LDT_verify
    use LDT_paramOptCheckMod, only: LDT_noahparmsOptChecks, &
                       LDT_gridOptChecks,LDT_soilsOptChecks
!
! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! the noahParms fraction datasets 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[noahParmssetup](\ref{noahParmssetup}) \newline
!    calls the registry to invoke the noahParms setup methods. 
!  \end{description}
!
!EOP
   implicit none
   integer, intent (in)      :: flag
   integer  :: n,i,c,r,m
   integer  :: rc
   real     :: temp
   logical  :: file_exists
   logical  :: check_data
   type(LDT_fillopts) :: tbot
   type(LDT_fillopts) :: slopetype
   real, allocatable  :: force_elev(:,:)
   character*50       :: tbot_proj
   character*50       :: slopetype_proj
   character*50       :: MMF_proj
   type(MMF_BCsReader):: MBR_FDEPTH, MBR_RECH, MBR_RIVERBED, MBR_WTD, MBR_HGT
   logical            :: run_mmf = .false.

! _____________________________________________________________________

   allocate( Noah_struc(LDT_rc%nnest) )
   do n=1,LDT_rc%nnest
      
      if(flag == 0) then
         call set_param_attribs(Noah_struc(n)%tbot, "TBOT",&
              units="K", &
              full_name="Noah LSM bottom temperature")
      elseif (flag == 1) then
         call set_param_attribs(Noah_struc(n)%tbot, "SOILTEMP",&
              units="K", &
              full_name="Noah LSM bottom temperature")
      endif
      
      call set_param_attribs(Noah_struc(n)%slopetype,"SLOPETYPE",&
            units="-", &
            full_name="Noah LSM slope type")

      if ((LDT_rc%lsm.eq."Noah-MP.3.6").or.                        &
          (LDT_rc%lsm.eq."Noah-MP.4.0.1")) then
         call set_param_attribs(Noah_struc(n)%pblh,"NOAHMP36_PBLH",&
               units="m", &
               full_name="Noah-MP LSM planetary boundary height")

         if (LDT_rc%lsm.eq."Noah-MP.4.0.1") then
            call ESMF_ConfigGetAttribute(LDT_config, run_mmf,label='Process Noah-MP-4.0.1 MMF groundwater parameters:', DEFAULT= .false., rc=rc)
            if (run_mmf) then
               call set_param_attribs(Noah_struc(n)%fdepth,"MMF_FDEPTH",&
                    units="m", &
                    full_name="transmissivity e-folding depth")
               call set_param_attribs(Noah_struc(n)%rechclim,"MMF_RECHCLIM",&
                    units="mm", &
                    full_name="climatological recharge")
               call set_param_attribs(Noah_struc(n)%riverbed,"MMF_RIVERBED",&
                    units="m", &
                    full_name="riverbed elevation")
               call set_param_attribs(Noah_struc(n)%eqwtd,"MMF_EQWTD",&
                    units="m", &
                    full_name="equilibrium water table depth")
               call set_param_attribs(Noah_struc(n)%areaxy,"AREAXY",&
                    units="km^2", &
                    full_name="area of the grid cell")
               call set_param_attribs(Noah_struc(n)%hgtm,"MMF_HGTM",&
                    units="m MSL", &
                    full_name="GMTED2010 30-arc-second topography height")
            endif
         endif
      endif
   enddo

! -- Slope type: --

   check_data = .false. 

   call ESMF_ConfigFindLabel(LDT_config,"Slope type data source:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,Noah_struc(n)%slopetype%source,rc=rc)
      call LDT_verify(rc,"Slope type data source: not defined")
      
      if( Noah_struc(n)%slopetype%source.eq."none" ) then 
         Noah_struc(n)%slopetype%selectOpt = 0
      endif
      if( Noah_struc(n)%slopetype%selectOpt.eq.1 ) then
         check_data = .true. 
         allocate(Noah_struc(n)%slopetype%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              Noah_struc(n)%slopetype%num_bins))
      endif
   enddo

   if(check_data) then 
        write(LDT_logunit,*)" - - - - - - - - - Slope Type Parameter - - - - - - - - - - - -"

    ! Read in slope type file config entries:
      call ESMF_ConfigFindLabel(LDT_config,"Slope type map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,Noah_struc(n)%slopetypefile,rc=rc)
      enddo
      call ESMF_ConfigGetAttribute(LDT_config,slopetype_proj,&
           label="Slope type map projection:",rc=rc)
      call LDT_verify(rc,'Slope type map projection: option not specified in the config file')
      Noah_struc(:)%slopetype_proj = slopetype_proj
      
      call ESMF_ConfigFindLabel(LDT_config,"Slope type spatial transform:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,Noah_struc(n)%slopetype_gridtransform,&
              rc=rc)
         call LDT_verify(rc,'Slope type spatial transform: option not specified in the config file')
      enddo
      
    ! Read in Slope type "fill" options:
      slopetype%filltype = "none"
      call ESMF_ConfigGetAttribute(LDT_config, slopetype%filltype, &
           label="Slope type fill option:",rc=rc)
      call LDT_verify(rc,"Slope type fill option: option not specified in the config file")
      
      if( slopetype%filltype == "neighbor" ) then
         call ESMF_ConfigGetAttribute(LDT_config, slopetype%fillradius, &
              label="Slope type fill radius:",rc=rc)
         call LDT_verify(rc,"Slope type fill radius: option not specified in the config file")
         
         call ESMF_ConfigGetAttribute(LDT_config, slopetype%fillvalue, &
              label="Slope type fill value:",rc=rc)
         call LDT_verify(rc,"Slope type fill value: option not specified in the config file")
         if( slopetype%fillvalue > 9. ) then
            slopetype%fillvalue = 9.  ! Set slopetype upper limit to 9 
         end if
      elseif( slopetype%filltype == "none" ) then
         write(LDT_logunit,*) " -- 'NONE' Parameter-Mask Agreement Option Selected for Slope type"
      else
         write(LDT_logunit,*) "[ERR] Fill option for Slope type is not valid: ",trim(slopetype%filltype)
         write(LDT_logunit,*) "  Please select one of these:  none or neighbor "
         write(LDT_logunit,*) "  Programming stopping ..."
         call LDT_endrun
      end if
      
      do n=1,LDT_rc%nnest
         if( index(Noah_struc(n)%slopetype%source,"Native").eq.0 .and. &
             index(Noah_struc(n)%slopetype%source,"CONSTANT").eq.0 ) then
            call LDT_readDomainConfigSpecs("Slope type", &
                     slopetype_proj, Noah_struc(n)%slopetype_gridDesc)
            if( slopetype_proj == "latlon" ) then
               call LDT_gridOptChecks( n, "Slope type", Noah_struc(n)%slopetype_gridtransform,&
                    slopetype_proj, Noah_struc(n)%slopetype_gridDesc(9) )
            endif
         endif
         call LDT_soilsOptChecks(n, "Slope type", &
                                 slopetype_proj, &
                                 Noah_struc(n)%slopetype_gridtransform )
      
      !- Read in slope type map (mostly used for Noah LSM):
         select case ( Noah_struc(n)%slopetype%source )
          case ( "NCEP_LIS" )
            call read_NCEP_slopetype(n,&
                      Noah_struc(n)%slopetype%value)
          case ( "NCEP_GFS" )
            call read_GFS_slopetype(n,&
                      Noah_struc(n)%slopetype%value)
          case ( "NCEP_Native" )
            call read_NCEPNative_slopetype(n,&
                      Noah_struc(n)%slopetype%value)
          case ( "CONSTANT" )
            call read_CONSTANT_slopetype(n,&
                      Noah_struc(n)%slopetype%value)
          case default
            write(LDT_logunit,*) "[WARN] Slopetype data source has not been selected."
            write(LDT_logunit,*) "  Your Noah LSM will not run without this parameter set."
            write(LDT_logunit,*) "  Please select one of the following: " 
            write(LDT_logunit,*) " -- NCEP_LIS, NCEP_Native, NCEP_GFS, CONSTANT "
            write(LDT_logunit,*) "Program stopping ..."
            call LDT_endrun 
         end select

       ! Fill where parameter values are missing compared to land/water mask:
         if( slopetype%filltype == "neighbor" ) then
            write(LDT_logunit,*) "Checking/filling mask values for: ", &
                 trim(Noah_struc(n)%slopetype%short_name)
            write(fill_logunit,*) "Checking/filling mask values for: ", &
                 trim(Noah_struc(n)%slopetype%short_name)
            slopetype%watervalue = 0.
            call LDT_discreteParam_Fill( n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
                 Noah_struc(n)%slopetype_gridtransform,                   &
                 Noah_struc(n)%slopetype%num_bins,                    &
                 Noah_struc(n)%slopetype%value, slopetype%watervalue, &
                 LDT_LSMparam_struc(n)%landmask2%value,               &
                 slopetype%filltype, slopetype%fillvalue, slopetype%fillradius )
         endif
     enddo

   end if  ! Slopetype selection check


!-- Bottom soil temperature (K) field:

   call ESMF_ConfigFindLabel(LDT_config,"Bottom temperature data source:",rc=rc)    
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,Noah_struc(n)%tbot%source,rc=rc)
      call LDT_verify(rc,'Bottom temperature data source: not specified')
   enddo

   check_data = .false. 
   do n=1,LDT_rc%nnest
      if(Noah_struc(n)%tbot%source.eq."none") then
         Noah_struc(n)%tbot%selectOpt = 0
      endif
      if(Noah_struc(n)%tbot%selectOpt.eq.1) then 
         check_data = .true. 
         allocate(Noah_struc(n)%tbot%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              Noah_struc(n)%tbot%vlevels))
      endif
   enddo

   if( check_data ) then 
       write(LDT_logunit,*)" - - - - - - - - - Bottom Temperature Parameter - - - - - - - - - - - -"

      call ESMF_ConfigFindLabel(LDT_config,"Bottom temperature map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,Noah_struc(n)%tbotfile,rc=rc)
         call LDT_verify(rc,'Bottom temperature map: not specified')
      enddo
      call ESMF_ConfigGetAttribute(LDT_config, tbot_proj,&
           label="Bottom temperature map projection:",rc=rc)
      call LDT_verify(rc,'Bottom temperature map projection: option not specified in the config file')
      Noah_struc(:)%tbot_proj = tbot_proj

      call ESMF_ConfigFindLabel(LDT_config,"Bottom temperature spatial transform:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,Noah_struc(n)%tbot_gridtransform,&
              rc=rc)
         call LDT_verify(rc,'Bottom temperature transform: option not specified in the config file')
      enddo

      Noah_struc(:)%tbot_topocorr = "none"
      call ESMF_ConfigFindLabel(LDT_config,"Bottom temperature topographic downscaling:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,Noah_struc(n)%tbot_topocorr,rc=rc)
         call LDT_verify(rc,'Bottom temperature topographic downscaling: not specified')
       ! Allow for mis-entered lapse-rate option entry:
         if( Noah_struc(n)%tbot_topocorr == "lapse_rate" ) Noah_struc(n)%tbot_topocorr="lapse-rate"
         if( Noah_struc(n)%tbot_topocorr == "lapse rate" ) Noah_struc(n)%tbot_topocorr="lapse-rate"
         if( Noah_struc(n)%tbot_topocorr == "Lapse-rate" ) Noah_struc(n)%tbot_topocorr="lapse-rate"
      enddo

    ! Read in "fill" option entries:
      tbot%filltype = "none"
      call ESMF_ConfigGetAttribute(LDT_config, tbot%filltype, &
           label="Bottom temperature fill option:",rc=rc)
      call LDT_verify(rc,"Bottom temperature fill option: option not specified in the config file")

      if( tbot%filltype == "neighbor" .or. tbot%filltype == "average" ) then
         call ESMF_ConfigGetAttribute(LDT_config, tbot%fillvalue, &
              label="Bottom temperature fill value:",rc=rc)
         call LDT_verify(rc,"Bottom temperature fill value: option not specified in the config file")

         call ESMF_ConfigGetAttribute(LDT_config, tbot%fillradius, &
              label="Bottom temperature fill radius:",rc=rc)
         call LDT_verify(rc,"Bottom temperature fill radius: option not specified in the config file")
      elseif( tbot%filltype == "none" ) then
         write(LDT_logunit,*) " -- 'NONE' Parameter-Mask Agreement Option Selected for Bottom Temperature"
      else
         write(LDT_logunit,*) "[ERR] Fill option for Bottom Temperature is not valid: ",trim(tbot%filltype)
         write(LDT_logunit,*) "  Please select one of these:  none, neighbor or average "
         write(LDT_logunit,*) "  Programming stopping ..."
         call LDT_endrun

      end if

    ! Don't need to read in "Native" grid extents/resolution, just for LIS inputs
      do n=1,LDT_rc%nnest
        if( index(Noah_struc(n)%tbot%source,"Native").eq.0  .and. &
            index(Noah_struc(n)%tbot%source,"ISLSCP1").eq.0 .and. &
            index(Noah_struc(n)%tbot%source,"CONSTANT").eq.0) then
           call LDT_readDomainConfigSpecs("Bottom temperature", &
                    tbot_proj, Noah_struc(n)%tbot_gridDesc)
           if( tbot_proj == "latlon" ) then
              call LDT_gridOptChecks( n, "Bottom temperature", Noah_struc(n)%tbot_gridtransform, &
                                      tbot_proj, Noah_struc(n)%tbot_gridDesc(9) )
           endif
        endif
        call LDT_noahparmsOptChecks( n, "Bottom temperature", tbot_proj, &
                                    Noah_struc(n)%tbot_gridtransform )
      enddo

      do n = 1, LDT_rc%nnest

       ! Read in Tbot File:
         select case ( Noah_struc(n)%tbot%source )
          case( "NCEP_LIS" )
            call read_NCEP_tbot(&
                      n,Noah_struc(n)%tbot%value(:,:,1))
          case( "NCEP_GFS" )
            call read_NCEP_GFS_tbot(&
                      n,Noah_struc(n)%tbot%value(:,:,1))
          case( "ISLSCP1" )
            call read_ISLSCP1_tbot(&
                      n,Noah_struc(n)%tbot%value(:,:,1))
          case( "CONSTANT" )
            call read_CONSTANT_tbot(&
                      n,Noah_struc(n)%tbot%value(:,:,1))
          case default
            write(LDT_logunit,*) "[WARN] Bottom temperature data source not selected."
            write(LDT_logunit,*) "  Your Noah LSM will not run without this parameter set."
            write(LDT_logunit,*) "  Please select one of the following: "
            write(LDT_logunit,*) " -- NCEP_LIS, ISLSCP1, NCEP_GFS, CONSTANT "
            write(LDT_logunit,*) "Program stopping ..."
            call LDT_endrun
         end select

       ! Fill where parameter values are missing compared to land/water mask:
         if( tbot%filltype == "neighbor" .or. &
              tbot%filltype == "average" ) then
            write(LDT_logunit,*) "Checking/filling mask values for: ", &
                 trim(Noah_struc(n)%tbot%short_name)
            write(fill_logunit,*) "Checking/filling mask values for: ", &
                 trim(Noah_struc(n)%tbot%short_name)
            tbot%watervalue = LDT_rc%udef
            call LDT_contIndivParam_Fill( n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
                 Noah_struc(n)%tbot_gridtransform,                         &
                 Noah_struc(n)%tbot%num_bins,                &
                 Noah_struc(n)%tbot%value, tbot%watervalue,  &
                 LDT_LSMparam_struc(n)%landmask2%value,              &
                 tbot%filltype, tbot%fillvalue, tbot%fillradius )
         endif

      !- Modify final Tbot output with elevation correction:
         if( Noah_struc(n)%tbot_topocorr == "lapse-rate" ) then
            if( LDT_LSMparam_struc(n)%elevation%selectOpt == 1 ) then
               write(LDT_logunit,*) "Performing lapse-rate correction to Tbot output."
               allocate(force_elev(LDT_rc%lnc(n),LDT_rc%lnr(n)))
               force_elev = 0.
               do r = 1, LDT_rc%lnr(n)
                  do c = 1, LDT_rc%lnc(n)
                     if( Noah_struc(n)%tbot%value(c,r,1)/=LDT_rc%udef ) &
                          Noah_struc(n)%tbot%value(c,r,1) =   &
                          Noah_struc(n)%tbot%value(c,r,1)     &
                          + (-0.0065)*(LDT_LSMparam_struc(n)%elevation%value(c,r,1)  &
                          - force_elev(c,r))
                  end do
               end do
               deallocate(force_elev)
            elseif( LDT_LSMparam_struc(n)%elevation%selectOpt == 0 ) then
               write(LDT_logunit,*) "Cannot perform lapse-rate correction to Tbot output,"
               write(LDT_logunit,*) " since no elevation/terrain map option was selected. "
               write(LDT_logunit,*) " Stopping ... "
               call LDT_endrun
            endif
         endif

      enddo
    end if

  !== Other Noah LSM related parameters ==

! -- Noah-MP Planetary Boundary Layer Height: --

   check_data = .false.
   if ((LDT_rc%lsm.eq."Noah-MP.3.6").or.(LDT_rc%lsm.eq."Noah-MP.4.0.1")) then

!   if(check_data) &! then
     write(LDT_logunit,*)" - - - - - - - - - Noah-MP Parameters - - - - - - - - - - - -"

     call ESMF_ConfigFindLabel(LDT_config,"Noah-MP PBL Height Value:",rc=rc)
     do n=1,LDT_rc%nnest
!        call ESMF_ConfigGetAttribute(LDT_config,Noah_struc(n)%pblh%source,rc=rc)
        call ESMF_ConfigGetAttribute(LDT_config,Noah_struc(n)%pblh_value,rc=rc)
        call LDT_verify(rc,"Noah-MP PBL Height Value: not defined")

        Noah_struc(n)%pblh%selectOpt = 1

        allocate(Noah_struc(n)%pblh%value(&
                 LDT_rc%lnc(n),LDT_rc%lnr(n),&
                 Noah_struc(n)%pblh%num_bins))

        Noah_struc(n)%pblh%value = Noah_struc(n)%pblh_value

      enddo
   endif

   !== MMF parameters
   RUN_MMF_SCHME: if (run_mmf) then
      
      check_data = .false.
      
      write(LDT_logunit,*)" - - - - - - - MMF Groundwater Parameters - - - - - - - - - -"
      
      do n = 1,LDT_rc%nnest         
         call ESMF_ConfigGetAttribute(LDT_config, Noah_struc(n)%mmf_fdepth_dir, label='MMF transmissivity dir:', rc=rc)
         call LDT_verify(rc,"MMF transmissivity dir: not defined")
         call ESMF_ConfigGetAttribute(LDT_config, Noah_struc(n)%mmf_rechclim_dir, label='MMF climatological recharge dir:',       RC=RC)
         call LDT_verify(rc,"MMF climatological recharge dir: not defined")
         call ESMF_ConfigGetAttribute(LDT_config, Noah_struc(n)%mmf_riverbed_dir, label='MMF riverbed elevation dir:',            RC=RC)
         call LDT_verify(rc,"MMF riverbed elevation dir: not defined")
         call ESMF_ConfigGetAttribute(LDT_config, Noah_struc(n)%mmf_eqwtd_dir   , label='MMF equilibrium water table depth dir:', RC=RC)
         call LDT_verify(rc,"MMF equilibrium water table depth dir: not defined")
         call ESMF_ConfigGetAttribute(LDT_config, Noah_struc(n)%mmf_hgtm_dir    , label='MMF HGT_M dir:',                         RC=RC)
         call LDT_verify(rc,"MMF HGT_M dir: not defined")
         check_data = .true.
      end do
      
      if (check_data) then
         
         call ESMF_ConfigGetAttribute(LDT_config, MMF_proj, label='MMF map projection:', rc=rc)
         call LDT_verify(rc,"MMF map projection not defined.")
         
         call ESMF_ConfigGetAttribute(LDT_config, MBR_FDEPTH%gap_fill%filltype,  label='FDEPTH fill option:', rc=rc)
         call LDT_verify(rc,"FDEPTH fill option not defined.")         
         call ESMF_ConfigGetAttribute(LDT_config, MBR_FDEPTH%gap_fill%fillradius,label='FDEPTH fill radius:', rc=rc)
         call LDT_verify(rc,"FDEPTH fill radius not defined.")  
         call ESMF_ConfigGetAttribute(LDT_config, MBR_FDEPTH%gap_fill%fillvalue, label='FDEPTH fill value:' , rc=rc)
         call LDT_verify(rc,"FDEPTH fill value not defined." )  
         call ESMF_ConfigGetAttribute(LDT_config, MBR_FDEPTH%water_value,        label='FDEPTH water value:', DEFAULT= 100., rc=rc)
         call LDT_verify(rc,"Suitable FDEPTH parameter value in lakes/waterbodies not defined." )        
         
         call ESMF_ConfigGetAttribute(LDT_config, MBR_RECH%gap_fill%filltype,  label='RECHCLIM fill option:',rc=rc)
         call LDT_verify(rc,"RECHCLIM fill option not defined.")       
         call ESMF_ConfigGetAttribute(LDT_config, MBR_RECH%gap_fill%fillradius,label='RECHCLIM fill radius:',rc=rc)
         call LDT_verify(rc,"RECHCLIM fill radius not defined.") 
         call ESMF_ConfigGetAttribute(LDT_config, MBR_RECH%gap_fill%fillvalue, label='RECHCLIM fill value:' ,rc=rc)
         call LDT_verify(rc,"RECHCLIM fill value not defined." )
         call ESMF_ConfigGetAttribute(LDT_config, MBR_RECH%water_value,        label='RECHCLIM water value:', DEFAULT= 0., rc=rc)
         call LDT_verify(rc,"Suitable RECHCLIM parameter value in lakes/waterbodies not defined." )
         
         call ESMF_ConfigGetAttribute(LDT_config, MBR_RIVERBED%gap_fill%filltype,  label='RIVERBED fill option:', DEFAULT= "neighbor", rc=rc)
         call LDT_verify(rc,"RIVERBED fill option not defined.")     
         call ESMF_ConfigGetAttribute(LDT_config, MBR_RIVERBED%gap_fill%fillradius,label='RIVERBED fill radius:', DEFAULT= 0., rc=rc)
         call LDT_verify(rc,"RIVERBED fill radius not defined.")  
         call ESMF_ConfigGetAttribute(LDT_config, MBR_RIVERBED%gap_fill%fillvalue, label='RIVERBED fill value:' , DEFAULT= 10000., rc=rc)
         call LDT_verify(rc,"RIVERBED fill value not defined." )
         call ESMF_ConfigGetAttribute(LDT_config, MBR_RIVERBED%water_value,        label='RIVERBED water value:', DEFAULT= -100., rc=rc)
         call LDT_verify(rc,"Suitable RIVERBED parameter value in lakes/waterbodies not defined." )
         
         call ESMF_ConfigGetAttribute(LDT_config, MBR_WTD%gap_fill%filltype,  label='EQWTD fill option:',rc=rc)
         call LDT_verify(rc,"EQWTD fill option not defined.")               
         call ESMF_ConfigGetAttribute(LDT_config, MBR_WTD%gap_fill%fillradius,label='EQWTD fill radius:',rc=rc)
         call LDT_verify(rc,"EQWTD fill radius not defined.")  
         call ESMF_ConfigGetAttribute(LDT_config, MBR_WTD%gap_fill%fillvalue, label='EQWTD fill value:' ,rc=rc)
         call LDT_verify(rc,"EQWTD fill value not defined." )
         call ESMF_ConfigGetAttribute(LDT_config, MBR_WTD%water_value,        label='EQWTD water value:', DEFAULT= 0., rc=rc)
         call LDT_verify(rc,"Suitable EQWTD parameter value in lakes/waterbodies not defined." )
         
         call ESMF_ConfigGetAttribute(LDT_config, MBR_HGT%gap_fill%filltype,  label='HGT_M fill option:',rc=rc)
         call LDT_verify(rc,"HGT_M fill option not defined.")               
         call ESMF_ConfigGetAttribute(LDT_config, MBR_HGT%gap_fill%fillradius,label='HGT_M fill radius:',rc=rc)
         call LDT_verify(rc,"HGT_M fill radius not defined.")  
         call ESMF_ConfigGetAttribute(LDT_config, MBR_HGT%gap_fill%fillvalue, label='HGT_M fill value:' ,rc=rc)
         call LDT_verify(rc,"HGT_M fill value not defined." )
          
         do n = 1,LDT_rc%nnest
            
            allocate (Noah_struc(n)%fdepth%value   (LDT_rc%lnc(n),LDT_rc%lnr(n),Noah_struc(n)%fdepth%num_bins  ))
            allocate (Noah_struc(n)%rechclim%value (LDT_rc%lnc(n),LDT_rc%lnr(n),Noah_struc(n)%rechclim%num_bins))
            allocate (Noah_struc(n)%riverbed%value (LDT_rc%lnc(n),LDT_rc%lnr(n),Noah_struc(n)%riverbed%num_bins))
            allocate (Noah_struc(n)%eqwtd%value    (LDT_rc%lnc(n),LDT_rc%lnr(n),Noah_struc(n)%eqwtd%num_bins   ))
            allocate (Noah_struc(n)%areaxy%value   (LDT_rc%lnc(n),LDT_rc%lnr(n),Noah_struc(n)%eqwtd%num_bins   ))
            allocate (Noah_struc(n)%hgtm%value     (LDT_rc%lnc(n),LDT_rc%lnr(n),Noah_struc(n)%hgtm%num_bins    ))
            
            call ESMF_ConfigGetAttribute(LDT_config, Noah_struc(n)%mmf_transform, label='MMF spatial transform:', rc=rc)
            call LDT_verify(rc,"MMF spatial transform method is not defined in the config file.")
            
            ! Read in index files in GEOGRID directories and create mapping
            
            call MBR_FDEPTH%mi     (n, MMF_proj, Noah_struc(n)%mmf_fdepth_dir)
            call MBR_RECH%mi       (n, MMF_proj, Noah_struc(n)%mmf_rechclim_dir, MBR_FDEPTH%MMF_mapping)
            call MBR_RIVERBED%mi   (n, MMF_proj, Noah_struc(n)%mmf_riverbed_dir, MBR_FDEPTH%MMF_mapping)
            call MBR_WTD%mi        (n, MMF_proj, Noah_struc(n)%mmf_eqwtd_dir   , MBR_FDEPTH%MMF_mapping)
            call MBR_HGT%mi        (n, MMF_proj, Noah_struc(n)%mmf_hgtm_dir    , MBR_FDEPTH%MMF_mapping)
            
            ! Read in variables fields
            
            call MBR_FDEPTH%mr   (n, trim(Noah_struc(n)%mmf_fdepth_dir  ), Noah_struc(n)%fdepth%value  , Noah_struc(n)%mmf_transform, Noah_struc(n)%fdepth%short_name  )
            call MBR_RECH%mr     (n, trim(Noah_struc(n)%mmf_rechclim_dir), Noah_struc(n)%rechclim%value, Noah_struc(n)%mmf_transform, Noah_struc(n)%rechclim%short_name)
            call MBR_RIVERBED%mr (n, trim(Noah_struc(n)%mmf_riverbed_dir), Noah_struc(n)%riverbed%value, Noah_struc(n)%mmf_transform, Noah_struc(n)%riverbed%short_name)
            call MBR_WTD%mr      (n, trim(Noah_struc(n)%mmf_eqwtd_dir   ), Noah_struc(n)%eqwtd%value   , Noah_struc(n)%mmf_transform, Noah_struc(n)%eqwtd%short_name   )
            call MBR_HGT%mr      (n, trim(Noah_struc(n)%mmf_hgtm_dir    ), Noah_struc(n)%hgtm%value    , Noah_struc(n)%mmf_transform, Noah_struc(n)%hgtm%short_name    )

            ! Ensure RIVERBED value in LAKES and WATER grid cells is same as MMF_HGTM
            where (Noah_struc(n)%riverbed%value == MBR_RIVERBED%water_value)
               Noah_struc(n)%riverbed%value = Noah_struc(n)%hgtm%value 
            endwhere

            ! Ensure RIVERBED value <= MMF_HGTM
            where (Noah_struc(n)%riverbed%value > Noah_struc(n)%hgtm%value)
               Noah_struc(n)%riverbed%value = Noah_struc(n)%hgtm%value 
            endwhere
            
            ! write areaXY            
            call cell_area (n,Noah_struc(n)%areaxy%value)
            
         end do       
      endif
   endif RUN_MMF_SCHME
   
 end subroutine NoahParms_init

 ! --------------------------------------------------------------------
 
 subroutine NoahParms_writeHeader(n,ftn,dimID)

    integer   :: n, rc 
    integer   :: ftn
    integer   :: dimID(3)
    logical   :: run_mmf = .false.
    
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
             Noah_struc(n)%tbot)
    
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
             Noah_struc(n)%slopetype)

    if ((LDT_rc%lsm.eq."Noah-MP.3.6").or.                        &
         (LDT_rc%lsm.eq."Noah-MP.4.0.1")) then
       call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
            Noah_struc(n)%pblh)
    endif
    if (LDT_rc%lsm.eq."Noah-MP.4.0.1") then
       call ESMF_ConfigGetAttribute(LDT_config, run_mmf,label='Process Noah-MP-4.0.1 MMF groundwater parameters:', DEFAULT= .false., rc=rc)
       if (run_mmf) then
          call LDT_writeNETCDFdataHeader(n,ftn,dimID,Noah_struc(n)%fdepth   ) 
          call LDT_writeNETCDFdataHeader(n,ftn,dimID,Noah_struc(n)%rechclim )
          call LDT_writeNETCDFdataHeader(n,ftn,dimID,Noah_struc(n)%riverbed )
          call LDT_writeNETCDFdataHeader(n,ftn,dimID,Noah_struc(n)%eqwtd    )
          call LDT_writeNETCDFdataHeader(n,ftn,dimID,Noah_struc(n)%areaxy   )
          call LDT_writeNETCDFdataHeader(n,ftn,dimID,Noah_struc(n)%hgtm     )     
       endif
    endif
  end subroutine NoahParms_writeHeader

  subroutine NoahParms_writeData(n,ftn)

    integer   :: n, rc 
    integer   :: ftn
    logical   :: run_mmf = .false.
    
    call LDT_writeNETCDFdata(n,ftn,Noah_struc(n)%tbot)

    call LDT_writeNETCDFdata(n,ftn,Noah_struc(n)%slopetype)

    if ((LDT_rc%lsm.eq."Noah-MP.3.6").or.                        &
        (LDT_rc%lsm.eq."Noah-MP.4.0.1")) then
        call LDT_writeNETCDFdata(n,ftn,Noah_struc(n)%pblh)
    endif
    if (LDT_rc%lsm.eq."Noah-MP.4.0.1") then
       call ESMF_ConfigGetAttribute(LDT_config, run_mmf,label='Process Noah-MP-4.0.1 MMF groundwater parameters:', DEFAULT= .false., rc=rc)
       if (run_mmf) then
          call LDT_writeNETCDFdata(n,ftn,Noah_struc(n)%fdepth   )
          call LDT_writeNETCDFdata(n,ftn,Noah_struc(n)%rechclim )
          call LDT_writeNETCDFdata(n,ftn,Noah_struc(n)%riverbed )
          call LDT_writeNETCDFdata(n,ftn,Noah_struc(n)%eqwtd    )
          call LDT_writeNETCDFdata(n,ftn,Noah_struc(n)%areaxy   )
          call LDT_writeNETCDFdata(n,ftn,Noah_struc(n)%hgtm     )    
       endif
    endif
  end subroutine NoahParms_writeData


!BOP
! !ROUTINE:  set_param_attribs
! \label{set_param_attribs}
!
! !INTERFACE:
  subroutine set_param_attribs(paramEntry, short_name, &
                               units, full_name )

! !DESCRIPTION:
!   This routine reads over the parameter attribute entries
!   in the param_attribs.txt file.
!
! !USES:
   type(LDT_paramEntry),intent(inout) :: paramEntry
   character(len=*),    intent(in)    :: short_name
   character(len=*),     optional     :: units
   character(len=*),     optional     :: full_name

   character(20) :: unit_temp
   character(100):: name_temp
! ____________________________________________________
    
   if(present(units)) then
      unit_temp = units
   else
      unit_temp = "none"
   endif
   if(present(full_name)) then
      name_temp = full_name
   else
      name_temp = trim(short_name)
   endif
   
   paramEntry%short_name = trim(short_name)
   paramEntry%vlevels = 1
   paramEntry%selectOpt = 1
   paramEntry%source = "Noah"
   paramEntry%units = unit_temp
   paramEntry%num_times = 1
   paramEntry%num_bins = 1
   paramEntry%standard_name = trim(name_temp)

  end subroutine set_param_attribs

end module Noah_parmsMod
