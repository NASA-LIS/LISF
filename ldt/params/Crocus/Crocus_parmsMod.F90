!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module Crocus_parmsMod
!BOP
!
! !MODULE: Crocus_parmsMod
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
!  23 Sep 2019: Mahdi Navari; This code is based on the subroutine 
!               lsm_paramsMod initially implemented by Sujay Kumar and 
!               modified by K. Arsenault. 
!  26 Jan 2021: Mahdi Navari; Modified for Crocus implementation
!
  use ESMF
  use LDT_coreMod
  use LDT_historyMod
  use LDT_paramDataMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_paramMaskCheckMod
  use LDT_glacierMod ! MN added for glacier fraction
  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: CrocusParms_init    !allocates memory for required structures
  public :: CrocusParms_writeHeader
  public :: CrocusParms_writeData

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: Crocus_struc

  type, public :: Crocus_type_dec

     real           :: tbot_gridDesc(20)
     character(len=LDT_CONST_PATH_LEN)  :: tbotfile
     character*50   :: tbot_gridtransform
     character*50   :: tbot_topocorr
     character(len=LDT_CONST_PATH_LEN)  :: glacierfracfile
     character*50   :: glacierfrac_gridtransform
     character*50   :: tbot_proj
     character*50   :: glacierfrac_proj


     ! -  Crocus SM-specific:
     type(LDT_paramEntry) :: tbot        ! Bottom temperature (Crocus)
     type(LDT_paramEntry) :: glacierfrac
 
  end type Crocus_type_dec

  type(Crocus_type_dec), allocatable :: Crocus_struc(:)

contains

!BOP
! 
! !ROUTINE: CrocusParms_init
! \label{CrocusParms_init}
! 
! !INTERFACE:
  subroutine CrocusParms_init
! !USES:
    use LDT_fileIOMod, only : LDT_readDomainConfigSpecs
    use LDT_logMod,    only : LDT_verify
    use LDT_paramOptCheckMod, only: LDT_noahparmsOptChecks, &
                       LDT_gridOptChecks,LDT_soilsOptChecks
!
! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! the CrocusParms fraction datasets 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[CrocusParmssetup](\ref{CrocusParmssetup}) \newline
!    calls the registry to invoke the CrocusParms setup methods. 
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
   real, allocatable  :: force_elev(:,:)
   character*50       :: tbot_proj
   character*50       :: glacierfrac_proj
   real, allocatable  :: localmask(:,:)
! _____________________________________________________________________

   allocate( Crocus_struc(LDT_rc%nnest) )
   do n=1,LDT_rc%nnest
      call set_param_attribs(Crocus_struc(n)%tbot, "TBOT",&
            units="K", &
            full_name="Crocus LSM bottom temperature")

      call set_param_attribs(Crocus_struc(n)%glacierfrac,"GlacierFraction",&
            units="-", &
            full_name="Crocus LSM Glacier Fraction")
   enddo



!-- Bottom soil temperature (K) field:

   call ESMF_ConfigFindLabel(LDT_config,"Bottom temperature data source:",rc=rc)    
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,Crocus_struc(n)%tbot%source,rc=rc)
      call LDT_verify(rc,'Bottom temperature data source: not specified')
   enddo

   check_data = .false. 
   do n=1,LDT_rc%nnest
      if(Crocus_struc(n)%tbot%source.eq."none") then
         Crocus_struc(n)%tbot%selectOpt = 0
      endif
      if(Crocus_struc(n)%tbot%selectOpt.eq.1) then 
         check_data = .true. 
         allocate(Crocus_struc(n)%tbot%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              Crocus_struc(n)%tbot%vlevels))
      endif
   enddo

   if( check_data ) then 
       write(LDT_logunit,*)" - - - - - - - - - Bottom Temperature Parameter - - - - - - - - - - - -"

      call ESMF_ConfigFindLabel(LDT_config,"Bottom temperature map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,Crocus_struc(n)%tbotfile,rc=rc)
         call LDT_verify(rc,'Bottom temperature map: not specified')
      enddo
      call ESMF_ConfigGetAttribute(LDT_config, tbot_proj,&
           label="Bottom temperature map projection:",rc=rc)
      call LDT_verify(rc,'Bottom temperature map projection: option not specified in the config file')
      Crocus_struc(:)%tbot_proj = tbot_proj

      call ESMF_ConfigFindLabel(LDT_config,"Bottom temperature spatial transform:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,Crocus_struc(n)%tbot_gridtransform,&
              rc=rc)
         call LDT_verify(rc,'Bottom temperature spatial transform: option not specified in the config file')
      enddo

      Crocus_struc(:)%tbot_topocorr = "none"
      call ESMF_ConfigFindLabel(LDT_config,"Bottom temperature topographic downscaling:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,Crocus_struc(n)%tbot_topocorr,rc=rc)
         call LDT_verify(rc,'Bottom temperature topographic downscaling: not specified')
       ! Allow for mis-entered lapse-rate option entry:
         if( Crocus_struc(n)%tbot_topocorr == "lapse_rate" ) Crocus_struc(n)%tbot_topocorr="lapse-rate"
         if( Crocus_struc(n)%tbot_topocorr == "lapse rate" ) Crocus_struc(n)%tbot_topocorr="lapse-rate"
         if( Crocus_struc(n)%tbot_topocorr == "Lapse-rate" ) Crocus_struc(n)%tbot_topocorr="lapse-rate"
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
        if( index(Crocus_struc(n)%tbot%source,"Native").eq.0  .and. &
            index(Crocus_struc(n)%tbot%source,"ISLSCP1").eq.0 .and. &
            index(Crocus_struc(n)%tbot%source,"CONSTANT").eq.0) then
           call LDT_readDomainConfigSpecs("Bottom temperature", &
                    tbot_proj, Crocus_struc(n)%tbot_gridDesc)
           if( tbot_proj == "latlon" ) then
              call LDT_gridOptChecks( n, "Bottom temperature", Crocus_struc(n)%tbot_gridtransform, &
                                      tbot_proj, Crocus_struc(n)%tbot_gridDesc(9) )
           endif
        endif
!        call LDT_noahparmsOptChecks( n, "Bottom temperature", tbot_proj, &
!                                    Noah_struc(n)%tbot_gridtransform )
      enddo

      do n = 1, LDT_rc%nnest

       ! Read in Tbot File:
         select case ( Crocus_struc(n)%tbot%source )
          case( "NCEP_LIS" )
            call read_NCEP_cro_tbot(&
                      n,Crocus_struc(n)%tbot%value(:,:,1))
          case( "NCEP_GFS" )
            call read_NCEP_GFS_cro_tbot(&
                      n,Crocus_struc(n)%tbot%value(:,:,1))
          case( "ISLSCP1" )
            call read_ISLSCP1_cro_tbot(&
                      n,Crocus_struc(n)%tbot%value(:,:,1))
          case( "CONSTANT" )
            call read_CONSTANT_tbot_cro(&
                      n,Crocus_struc(n)%tbot%value(:,:,1))
          case default
            write(LDT_logunit,*) "[WARN] Bottom temperature data source not selected."
            write(LDT_logunit,*) "  Your Crocus LSM will not run without this parameter set."
            write(LDT_logunit,*) "  Please select one of the following: "
            write(LDT_logunit,*) " -- NCEP_LIS, ISLSCP1, NCEP_GFS, CONSTANT "
            write(LDT_logunit,*) "Program stopping ..."
            call LDT_endrun
         end select

       ! Fill where parameter values are missing compared to land/water mask:
         if( tbot%filltype == "neighbor" .or. &
              tbot%filltype == "average" ) then
            write(LDT_logunit,*) "Checking/filling mask values for: ", &
                 trim(Crocus_struc(n)%tbot%short_name)
            write(fill_logunit,*) "Checking/filling mask values for: ", &
                 trim(Crocus_struc(n)%tbot%short_name)
            tbot%watervalue = LDT_rc%udef
            call LDT_contIndivParam_Fill( n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
                 Crocus_struc(n)%tbot_gridtransform,                         &
                 Crocus_struc(n)%tbot%num_bins,                &
                 Crocus_struc(n)%tbot%value, tbot%watervalue,  &
                 LDT_LSMparam_struc(n)%landmask2%value,              &
                 tbot%filltype, tbot%fillvalue, tbot%fillradius )
         endif

      !- Modify final Tbot output with elevation correction:
         if( Crocus_struc(n)%tbot_topocorr == "lapse-rate" ) then
            if( LDT_LSMparam_struc(n)%elevation%selectOpt == 1 ) then
               write(LDT_logunit,*) "Performing lapse-rate correction to Tbot output."
               allocate(force_elev(LDT_rc%lnc(n),LDT_rc%lnr(n)))
               force_elev = 0.
               do r = 1, LDT_rc%lnr(n)
                  do c = 1, LDT_rc%lnc(n)
                     if( Crocus_struc(n)%tbot%value(c,r,1)/=LDT_rc%udef ) &
                          Crocus_struc(n)%tbot%value(c,r,1) =   &
                          Crocus_struc(n)%tbot%value(c,r,1)     &
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
  end subroutine CrocusParms_init

  subroutine CrocusParms_writeHeader(n,ftn,dimID)

    integer   :: n 
    integer   :: ftn
    integer   :: dimID(3)

    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
             Crocus_struc(n)%tbot)
  end subroutine CrocusParms_writeHeader

  subroutine CrocusParms_writeData(n,ftn)

    integer   :: n 
    integer   :: ftn

    call LDT_writeNETCDFdata(n,ftn,Crocus_struc(n)%tbot)
  end subroutine CrocusParms_writeData


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
   paramEntry%source = "Crocus"
   paramEntry%units = unit_temp
   paramEntry%num_times = 1
   paramEntry%num_bins = 1
   paramEntry%standard_name = trim(name_temp)

  end subroutine set_param_attribs

end module Crocus_parmsMod
