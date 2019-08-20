!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
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
!
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
  public :: NoahParms_init    !allocates memory for required structures
  public :: NoahParms_writeHeader
  public :: NoahParms_writeData

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: Noah_struc

  type, public :: noah_type_dec

     real           :: tbot_gridDesc(20)
     character*140  :: tbotfile
     character*50   :: tbot_gridtransform
     character*50   :: tbot_topocorr

     real           :: slopetype_gridDesc(20)
     character*140  :: slopetypefile
     character*50   :: slopetype_gridtransform

     character*50   :: tbot_proj
     character*50   :: slopetype_proj

     real           :: pblh_value

     ! -  Noah LSM-specific:
     type(LDT_paramEntry) :: tbot        ! Bottom temperature (Noah)
     type(LDT_paramEntry) :: slopetype   ! Slope type index (Noah)
     type(LDT_paramEntry) :: pblh        ! Planetary Boundary Layer Height (Noah-MP)

  end type noah_type_dec

  type(noah_type_dec), allocatable :: Noah_struc(:)

contains

!BOP
! 
! !ROUTINE: NoahParms_init
! \label{NoahParms_init}
! 
! !INTERFACE:
  subroutine NoahParms_init
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

! _____________________________________________________________________

   allocate( Noah_struc(LDT_rc%nnest) )
   do n=1,LDT_rc%nnest
      call set_param_attribs(Noah_struc(n)%tbot, "TBOT",&
            units="K", &
            full_name="Noah LSM bottom temperature")
      call set_param_attribs(Noah_struc(n)%slopetype,"SLOPETYPE",&
            units="-", &
            full_name="Noah LSM slope type")

      if ((LDT_rc%lsm.eq."Noah-MP.3.6").or.                        &
          (LDT_rc%lsm.eq."Noah-MP.4.0.1")) then
         call set_param_attribs(Noah_struc(n)%pblh,"NOAHMP36_PBLH",&
               units="m", &
               full_name="Noah-MP LSM planetary boundary height")
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


  end subroutine NoahParms_init

  subroutine NoahParms_writeHeader(n,ftn,dimID)

    integer   :: n 
    integer   :: ftn
    integer   :: dimID(3)

    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
             Noah_struc(n)%tbot)
    
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
             Noah_struc(n)%slopetype)

    if ((LDT_rc%lsm.eq."Noah-MP.3.6").or.                        &
        (LDT_rc%lsm.eq."Noah-MP.4.0.1")) then
        call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
                 Noah_struc(n)%pblh)
    endif

  end subroutine NoahParms_writeHeader

  subroutine NoahParms_writeData(n,ftn)

    integer   :: n 
    integer   :: ftn

    call LDT_writeNETCDFdata(n,ftn,Noah_struc(n)%tbot)

    call LDT_writeNETCDFdata(n,ftn,Noah_struc(n)%slopetype)

    if ((LDT_rc%lsm.eq."Noah-MP.3.6").or.                        &
        (LDT_rc%lsm.eq."Noah-MP.4.0.1")) then
        call LDT_writeNETCDFdata(n,ftn,Noah_struc(n)%pblh)
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
