!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
module LDT_gfracMod
!BOP
!
! !MODULE: LDT_gfracMod
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
!
!$ use omp_lib
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
  public :: LDT_gfrac_readParamSpecs
  public :: LDT_greenness_init    !allocates memory for required structures
  public :: LDT_greenness_writeHeader
  public :: LDT_greenness_writeData

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LDT_gfrac_struc

  type, public :: gfrac_type_dec
     real                   :: gfrac_gridDesc(20)
     character*50           :: gfrac_proj
     character*50           :: gfrac_gridtransform
     character*100          :: gfracdir
     character*140          :: gfracfile
     character*20           :: gfracInterval

     character*140   :: shdmaxfile
     character*140   :: shdminfile

     type(LDT_paramEntry) :: gfrac       ! Clim.-based greenness fraction
     type(LDT_paramEntry) :: shdmin      ! Min. Greenness fraction (@pixel)
     type(LDT_paramEntry) :: shdmax      ! Max. Greenness fraction (@pixel)
  end type gfrac_type_dec

  type(gfrac_type_dec), allocatable :: LDT_gfrac_struc(:)


contains

  subroutine LDT_gfrac_readParamSpecs
    
    character*100     :: source
    integer           :: rc
    integer           :: n
    
    allocate(LDT_gfrac_struc(LDT_rc%nnest))

  ! Greenness fraction files:
    call ESMF_ConfigFindLabel(LDT_config,"Greenness data source:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
!       call LDT_warning(rc,"Greenness data source: not defined")
       call LDT_set_param_attribs(rc,LDT_gfrac_struc(n)%gfrac,&
            "GREENNESS",source)
       call LDT_set_param_attribs(rc,LDT_gfrac_struc(n)%shdmax,&
            "SHDMAX",source)
       call LDT_set_param_attribs(rc,LDT_gfrac_struc(n)%shdmin,&
            "SHDMIN",source)
    enddo
  ! LSM-required parameter check:
    if( index(LDT_rc%lsm,"Noah") == 1 .or. &
        index(LDT_rc%lsm,"CLSM") == 1 .or. &
        index(LDT_rc%lsm,"RDHM") == 1 .or. &
        index(LDT_rc%lsm,"SACHTET") == 1 ) then
      if( rc /= 0 ) then
         call LDT_warning(rc,"WARNING: Greenness data source: not defined")
      endif
    endif

  end subroutine LDT_gfrac_readParamSpecs

!BOP
! 
! !ROUTINE: LDT_greenness_init
! \label{LDT_greenness_init}
! 
! !INTERFACE:
  subroutine LDT_greenness_init
! !USES:
    use LDT_fileIOMod, only : LDT_readDomainConfigSpecs
    use LDT_logMod,    only : LDT_verify
    use LDT_paramOptCheckMod, only: LDT_gfracOptChecks, &
                       LDT_gridOptChecks

! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! the greenness fraction datasets 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[gfracsetup](\ref{gfracsetup}) \newline
!    calls the registry to invoke the gfrac setup methods. 
!  \end{description}
!
!EOP
   implicit none
   integer     :: n, k, c, r
   integer     :: rc
   character*3 :: months(12)
   character*3 :: sacmonths(12)
   data months /'jan','feb','mar','apr','may','jun','jul','aug',&
                'sep','oct','nov','dec'/
   data sacmonths /'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug',&
                   'Sep','Oct','Nov','Dec'/
   logical     :: file_exists
   type(ESMF_Config)  :: rdhmconsts_table
   type(LDT_fillopts) :: gfrac
   type(LDT_fillopts) :: shdmax
   type(LDT_fillopts) :: shdmin
   logical            :: gfrac_select, check_data1, check_data2
   logical            :: calc_minmaxgfrac

   real, allocatable            :: gfrac_gridDesc(:,:)
   character*50             :: gfrac_proj
   character*50,  allocatable   :: gfrac_gridtransform(:)
   character*100, allocatable   :: gfracdir(:)
   character*140, allocatable   :: gfracfile(:)
   character*20,  allocatable   :: gfracInterval(:)

! _____________________________________________________________

   gfrac_select = .false. 
   do n=1,LDT_rc%nnest
      if(LDT_gfrac_struc(n)%gfrac%selectOpt.gt.0) then 
         gfrac_select = .true. 
      endif
   enddo

   if(gfrac_select) then
      write(LDT_logunit,*)" - - - - - - - - - Greenness Fraction Parameters - - - - - - - - -"
   endif

   allocate(gfrac_gridDesc(LDT_rc%nnest,20))
   allocate(gfrac_gridtransform(LDT_rc%nnest))         
   allocate(gfracfile(LDT_rc%nnest))
   allocate(gfracdir(LDT_rc%nnest))
   allocate(gfracInterval(LDT_rc%nnest))

   do n=1,LDT_rc%nnest

      if(LDT_gfrac_struc(n)%gfrac%selectOpt.gt.0) then 
 
         call set_gfrac_attribs( n, LDT_gfrac_struc(n)%gfrac%source )
         LDT_gfrac_struc(n)%gfrac%vlevels  = LDT_gfrac_struc(n)%gfrac%num_times

         allocate(LDT_gfrac_struc(n)%gfrac%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              LDT_gfrac_struc(n)%gfrac%vlevels))       
      endif
      if(LDT_gfrac_struc(n)%shdmax%selectOpt.gt.0) then 
         LDT_gfrac_struc(n)%shdmax%vlevels = LDT_gfrac_struc(n)%shdmax%num_times
         allocate(LDT_gfrac_struc(n)%shdmax%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              LDT_gfrac_struc(n)%shdmax%vlevels))       
      endif
      if(LDT_gfrac_struc(n)%shdmin%selectOpt.gt.0) then 
         LDT_gfrac_struc(n)%shdmin%vlevels = LDT_gfrac_struc(n)%shdmin%num_times
         allocate(LDT_gfrac_struc(n)%shdmin%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              LDT_gfrac_struc(n)%shdmin%vlevels))       
      endif

   enddo

   !- Greenness fraction:
   if( gfrac_select ) then 

      call ESMF_ConfigFindLabel(LDT_config,"Greenness fraction map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,gfracdir(n),rc=rc)
         call LDT_verify(rc,'Greenness fraction map: not specified')
         LDT_gfrac_struc(n)%gfracdir = gfracdir(n)
      enddo

      call ESMF_ConfigFindLabel(LDT_config,"Greenness climatology interval:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,gfracInterval(n),rc=rc)
         call LDT_verify(rc,'Greenness climatology interval: not specified')
         if( trim(gfracInterval(n)) .ne. "monthly" ) then
            write(LDT_logunit,*) "ERR: 'monthly' Greenness fraction interval option specified."
            write(LDT_logunit,*) "    Set ... Greenness climatology interval:  monthly"
            write(LDT_logunit,*) "    Stopping."
            call LDT_endrun
         endif
         if( gfracInterval(n) == "monthly" .and. &
              LDT_gfrac_struc(n)%gfrac%vlevels .ne. 12 )then
            write(LDT_logunit,*) "ERR: The 'monthly' Greenness climatology interval option "
            write(LDT_logunit,*) "   should have '12' in the parameter attribs table."
            write(LDT_logunit,*) "   Please change to '12' there.  Stopping."
            call LDT_endrun
         endif
         LDT_gfrac_struc(n)%gfracInterval = gfracInterval(n)
      enddo

    ! Set units and full names:
      do n=1,LDT_rc%nnest
         LDT_gfrac_struc(n)%gfrac%units="-"
         call setGfracParmsFullnames( n, "gfrac", &
                 LDT_gfrac_struc(n)%gfrac%source )
      enddo

      call ESMF_ConfigGetAttribute(LDT_config,gfrac_proj,&
           label="Greenness map projection:",rc=rc)
      call LDT_verify(rc,'gfrac projection: option not specified in the config file')

      LDT_gfrac_struc(:)%gfrac_proj = gfrac_proj

      call ESMF_ConfigFindLabel(LDT_config,"Greenness spatial transform:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,gfrac_gridtransform(n),&
              rc=rc)
         call LDT_verify(rc,'Greenness spatial transform: option not specified in the config file')
         LDT_gfrac_struc(n)%gfrac_gridtransform = gfrac_gridtransform(n)
      enddo
   endif

!- Check if greenness min and max fraction fields are selected:
   check_data1 = .false. 
   check_data2 = .false.
   do n=1,LDT_rc%nnest
      if(LDT_gfrac_struc(n)%shdmax%selectOpt.gt.0 ) then 
         check_data1 = .true. 
      endif
      if(LDT_gfrac_struc(n)%shdmin%selectOpt.gt.0 ) then
         check_data2 = .true.
      endif
   enddo

!- Decide to either calculate or read-in min/max greenness values:
   calc_minmaxgfrac = .false.
   if( check_data1 .or. check_data2 ) then
     call ESMF_ConfigGetAttribute(LDT_config,calc_minmaxgfrac,&
               label="Calculate min-max greenness fraction:",rc=rc)
     call LDT_verify(rc,'Calculate min-max greenness fraction: not specified')
     if( .NOT. calc_minmaxgfrac ) then
       if( check_data1 ) then 
        call ESMF_ConfigFindLabel(LDT_config,"Greenness maximum map:",rc=rc)
        do n=1,LDT_rc%nnest
           call ESMF_ConfigGetAttribute(LDT_config,LDT_gfrac_struc(n)%shdmaxfile,rc=rc)
           call LDT_verify(rc,'Greenness maximum map: not specified')
        enddo
       endif
       if( check_data2 ) then 
        call ESMF_ConfigFindLabel(LDT_config,"Greenness minimum map:",rc=rc)
        do n=1,LDT_rc%nnest
           call ESMF_ConfigGetAttribute(LDT_config,LDT_gfrac_struc(n)%shdminfile,rc=rc)
           call LDT_verify(rc,'Greenness minimum map: not specified')
        enddo
       endif
     endif
   ! Set units and full names:
     do n=1,LDT_rc%nnest
        LDT_gfrac_struc(n)%shdmin%units="-"
        call setGfracParmsFullnames( n, "shdmin", &
                LDT_gfrac_struc(n)%shdmin%source )
        LDT_gfrac_struc(n)%shdmax%units="-"
        call setGfracParmsFullnames( n, "shdmax", &
                LDT_gfrac_struc(n)%shdmax%source )
     enddo
   endif

   !- Mask-parameter agreement options:
   if( gfrac_select .or. check_data1 .or. check_data2 ) then  
      gfrac%filltype = "none"
      call ESMF_ConfigGetAttribute(LDT_config, gfrac%filltype, &
           label="Greenness fill option:",rc=rc)
      call LDT_verify(rc,"Greenness fill option: option not specified in the config file")

      if( gfrac%filltype == "average" .or. gfrac%filltype == "neighbor" ) then
         call ESMF_ConfigGetAttribute(LDT_config, gfrac%fillradius, &
              label="Greenness fill radius:",rc=rc)
         call LDT_verify(rc,"Greenness fill radius: option not specified in the config file")

         call ESMF_ConfigGetAttribute(LDT_config, gfrac%fillvalue, &
              label="Greenness fill value:",rc=rc)
         call LDT_verify(rc,"Greenness fill value: option not specified in the config file")

         if( LDT_gfrac_struc(1)%shdmax%selectOpt == 1 ) then
            call ESMF_ConfigGetAttribute(LDT_config, shdmax%fillvalue, &
                 label="Greenness maximum fill value:",rc=rc)
            call LDT_verify(rc,"Greenness maximum fill value: option not specified in the config file")
         endif
         if( LDT_gfrac_struc(1)%shdmin%selectOpt == 1 ) then
            call ESMF_ConfigGetAttribute(LDT_config, shdmin%fillvalue, &
                 label="Greenness minimum fill value:",rc=rc)
            call LDT_verify(rc,"Greenness minimum fill value: option not specified in the config file")
         endif
      elseif( gfrac%filltype == "none" ) then
         write(LDT_logunit,*) "[INFO] 'NONE' Parameter-Mask Agreement Option Selected for Greenness Fraction"
      else
         write(LDT_logunit,*) "[ERR] Fill option for Greenness fraction is not valid: ",trim(gfrac%filltype)
         write(LDT_logunit,*) "  Please select one of these:  none, neighbor or average "
         write(LDT_logunit,*) "  Programming stopping ..."
         call LDT_endrun
      end if

   endif

   !- Greenness Fraction:
   do n = 1, LDT_rc%nnest

      if( LDT_gfrac_struc(n)%gfrac%selectOpt .gt. 0) then 
       ! Don't need to read in "Native" grid extents/resolution, just for LIS inputs
         if( index(LDT_gfrac_struc(n)%gfrac%source,"Native").eq.0  .and. &
             index(LDT_gfrac_struc(n)%gfrac%source,"CONSTANT").eq.0  ) then
           call LDT_readDomainConfigSpecs("Greenness", gfrac_proj, gfrac_gridDesc)
           if( gfrac_proj == "latlon" ) then
             call LDT_gridOptChecks( n, "Greenness", gfrac_gridtransform(n), &
                                     gfrac_proj, gfrac_gridDesc(n,9) )
           endif
        endif
        LDT_gfrac_struc(n)%gfrac_gridDesc = gfrac_gridDesc(n,:)
     endif
     
  enddo
  
   do n = 1, LDT_rc%nnest

      if( LDT_gfrac_struc(n)%gfrac%selectOpt .gt. 0) then 

        call LDT_gfracOptChecks( "Greenness", gfrac_proj, gfrac_gridtransform(n) )

      !- Monthly (interval) files:
         if( gfracInterval(n) == "monthly" ) then 
             LDT_rc%monthlyData(n) = .true.

         !- Read single-file monthly clim greenness fraction: 
            if( trim(LDT_gfrac_struc(n)%gfrac%source) == "CLSMF2.5" ) then

               gfracfile(n) = trim(gfracdir(n))
               LDT_gfrac_struc(n)%gfracfile = gfracfile(n)

               write(LDT_logunit,*) "Reading single-file, monthly climatologies for: "&
                    //trim(gfracfile(n))
               call readgfrac( trim(LDT_gfrac_struc(n)%gfrac%source)//char(0),&
                    n, LDT_gfrac_struc(n)%gfrac%value, &
                    LDT_LSMparam_struc(n)%landmask%value )
               write(LDT_logunit,*) "Done reading file - "//trim(gfracfile(n))

         !- Read multi-file monthly clim greenness fraction: 
            else
               do k = 1, LDT_gfrac_struc(n)%gfrac%vlevels

                  if( trim(LDT_gfrac_struc(n)%gfrac%source) == "SACHTET.3.5.6" ) then
                     gfracfile(n) = trim(gfracdir(n))//"_"//&
                          trim(sacmonths(k))//".gz"
                  elseif( trim(LDT_gfrac_struc(n)%gfrac%source).eq."NCEP_LIS" ) then
                     gfracfile(n) = trim(gfracdir(n))//"."//&
                          trim(months(k))//".1gd4r"
                  elseif( trim(LDT_gfrac_struc(n)%gfrac%source).eq."NCEP_Native" )then
                     gfracfile(n) = trim(gfracdir(n))//"_"//&
                          trim(months(k))//".asc"
                  elseif( trim(LDT_gfrac_struc(n)%gfrac%source).eq."CONSTANT" )then
                     gfracfile(n) = trim(gfracdir(n))
                  endif
                  LDT_gfrac_struc(n)%gfracfile = gfracfile(n)

                  write(LDT_logunit,*) "Reading "//trim(gfracfile(n))
                  call readgfrac( trim(LDT_gfrac_struc(n)%gfrac%source)//char(0),&
                       n, LDT_gfrac_struc(n)%gfrac%value(:,:,k) )
                  write(LDT_logunit,*) "Done reading "//trim(gfracfile(n))

               enddo
            end if
         endif   ! End monthly interval 

         if( gfrac%filltype == "average" .or. gfrac%filltype == "neighbor" ) then
            write(LDT_logunit,*) "Checking/filling mask values for: ", &
                 trim(LDT_gfrac_struc(n)%gfrac%short_name)
            write(fill_logunit,*) "Checking/filling mask values for: ", &
                 trim(LDT_gfrac_struc(n)%gfrac%short_name)
            gfrac%watervalue = LDT_rc%udef
            call LDT_contIndivParam_Fill(n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
                 gfrac_gridtransform(n), LDT_gfrac_struc(n)%gfrac%num_times, &
                 LDT_gfrac_struc(n)%gfrac%value, gfrac%watervalue,     &
                 LDT_LSMparam_struc(n)%landmask2%value, gfrac%filltype,   &
                 gfrac%fillvalue, gfrac%fillradius )
         endif
      endif   ! End Gfrac map read

   !- Max greenness fraction:
      if( LDT_gfrac_struc(n)%shdmax%selectOpt == 1 ) then 
 
     !- Calculate maximum greenness values:
        if( calc_minmaxgfrac ) then
          do r = 1, LDT_rc%lnr(n)
            do c = 1, LDT_rc%lnc(n)
              LDT_gfrac_struc(n)%shdmax%value(c,r,1) = &
                  maxval(LDT_gfrac_struc(n)%gfrac%value(c,r,:))
            enddo
          enddo
     !- Read maximum greenness fraction file:
        else
          call readshdmax(trim(LDT_gfrac_struc(n)%gfrac%source)//char(0),&
                          n,LDT_gfrac_struc(n)%shdmax%value(:,:,1))
        endif
     !- Perform mask-parameter consistency "fill" options:
        if( gfrac%filltype == "average" .or. gfrac%filltype == "neighbor" ) then
           write(LDT_logunit,*) "Checking/filling mask values for: ", &
                trim(LDT_gfrac_struc(n)%shdmax%short_name)
           write(fill_logunit,*) "Checking/filling mask values for: ", &
                trim(LDT_gfrac_struc(n)%shdmax%short_name)
           gfrac%watervalue = LDT_rc%udef
           call LDT_contIndivParam_Fill(n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
                gfrac_gridtransform(n), LDT_gfrac_struc(n)%shdmax%num_times, &
                LDT_gfrac_struc(n)%shdmax%value, gfrac%watervalue,      &
                LDT_LSMparam_struc(n)%landmask2%value, gfrac%filltype,      &
                shdmax%fillvalue, gfrac%fillradius )
        endif
      endif

   !- Min greenness fraction:
      if( LDT_gfrac_struc(n)%shdmin%selectOpt == 1 ) then 

     !- Calculate minimum greenness values:
        if( calc_minmaxgfrac ) then
          do r = 1, LDT_rc%lnr(n)
            do c = 1, LDT_rc%lnc(n)
               LDT_gfrac_struc(n)%shdmin%value(c,r,1) = &
                   minval(LDT_gfrac_struc(n)%gfrac%value(c,r,:))
            enddo
          enddo
     !- Read minimum greenness fraction file:
        else
       !- Read minimum greenness fraction file:
          call readshdmin(trim(LDT_gfrac_struc(n)%gfrac%source)//char(0),&
                          n,LDT_gfrac_struc(n)%shdmin%value(:,:,1))
        endif
     !- Perform mask-parameter consistency "fill" options:
        if( gfrac%filltype == "average" .or. gfrac%filltype == "neighbor" ) then
           write(LDT_logunit,*) "Checking/filling mask values for: ", &
                trim(LDT_gfrac_struc(n)%shdmin%short_name)
           write(fill_logunit,*) "Checking/filling mask values for: ", &
                trim(LDT_gfrac_struc(n)%shdmin%short_name)
           gfrac%watervalue = LDT_rc%udef
           call LDT_contIndivParam_Fill(n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
                gfrac_gridtransform(n), LDT_gfrac_struc(n)%shdmin%num_times, &
                LDT_gfrac_struc(n)%shdmin%value, gfrac%watervalue,      &
                LDT_LSMparam_struc(n)%landmask2%value, gfrac%filltype,      &
                shdmin%fillvalue, gfrac%fillradius )
        endif
      endif

   enddo    ! End nest loop

  end subroutine LDT_greenness_init

  subroutine LDT_greenness_writeHeader(n,ftn,dimID,monthID)

    integer     :: n 
    integer     :: ftn
    integer     :: dimID(3)
    integer     :: monthID

    integer     :: t_dimID(3)

    if(LDT_gfrac_struc(n)%gfrac%selectOpt.gt.0) then 
       if(LDT_gfrac_struc(n)%gfracInterval.eq."monthly") then !monthly
          t_dimID(1) = dimID(1)
          t_dimID(2) = dimID(2)
          t_dimID(3) = monthID       
   
          call LDT_writeNETCDFdataHeader(n,ftn,t_dimID,&
               LDT_gfrac_struc(n)%gfrac)
          
          call LDT_writeNETCDFdataHeader(n,ftn,t_dimID,&
               LDT_gfrac_struc(n)%shdmin)
          
          call LDT_writeNETCDFdataHeader(n,ftn,t_dimID,&
               LDT_gfrac_struc(n)%shdmax)
       endif

       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"GREENNESS_DATA_INTERVAL", &
            LDT_gfrac_struc(n)%gfracInterval))
    endif
  end subroutine LDT_greenness_writeHeader

  subroutine LDT_greenness_writeData(n,ftn)

    integer     :: n 
    integer     :: ftn

    if(LDT_gfrac_struc(n)%gfrac%selectOpt.gt.0) then 
       call LDT_writeNETCDFdata(n,ftn,LDT_gfrac_struc(n)%gfrac)
       call LDT_writeNETCDFdata(n,ftn,LDT_gfrac_struc(n)%shdmin)
       call LDT_writeNETCDFdata(n,ftn,LDT_gfrac_struc(n)%shdmax)
    endif
  end subroutine LDT_greenness_writeData

end module LDT_gfracMod
