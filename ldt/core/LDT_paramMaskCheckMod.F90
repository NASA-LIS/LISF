!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_paramMaskCheckMod
!BOP
!
! !MODULE: LDT_paramMaskCheckMod
!
! !DESCRIPTION:
!  The following code provides means for checking and ensuring 
!   consistency between each selectec parameter and the read-in
!   or created binary land/water mask.  (Can be expanded later
!   to include surface type map).
!
!  \subsubsection{Overview}
!  This module contains routines that:
!  - perform parameter and mask consistency checks 
!  - ensures that each selected parameter has a valid value for 
!     corresponding mask-land values
!  - Performs nearest neighbor, average or dominant type for filling,
!     depending on the parameter type
!    (this section can be expanded with future releases)
! 
!  \begin{description}
!   \item[binary land/water mask]
!   \item[selected parameter]
!  \end{description}
!
! !REVISION HISTORY:
!
!  09 Aug 2012: KR Arsenault; Initial implementation
!
  use LDT_historyMod
  use LDT_paramDataMod
  use LDT_coreMod
  use LDT_logMod

  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LDT_ParamMaskFill_Log       ! Open Parameter-Mask Fill Log File
  public :: LDT_discreteParam_Fill      ! Fill routine for discrete parameters
  public :: LDT_contIndivParam_Fill     ! Fill for independent continuous parms
  public :: LDT_contTileParam_Fill      ! Fill for tiled continuous parms
!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  integer,public :: fill_logunit
  type, public :: LDT_fillopts
    character*50 :: filltype        ! Fill option type (none,neighbor,average)
    real         :: fillradius      ! Radius to search for neighboring values
    real         :: fillvalue       ! Value to fill in missing parameter value
    real         :: watervalue      ! Water or undefined value in parameter 
    real         :: fillvalue_antarctica ! EMK fillvalue south of 60S
    logical      :: force_exclude_water ! EMK Exclude water points during fill
  end type LDT_fillopts

!EOP

contains

!BOP
! 
! !ROUTINE: LDT_ParamMaskFill_Log
! \label{LDT_ParamMaskFill_Log}
! 
! !INTERFACE:
 subroutine LDT_ParamMaskFill_Log
!
! !DESCRIPTION:
! 
!  This routine opens a log file to write output
!   from parameter-mask fill options and points filled.
!
   fill_logunit = LDT_getNextUnitNumber()
!   open( fill_logunit, file = "MaskParamFill.log", form="formatted" )
   open( fill_logunit, file = LDT_rc%mpfillfile, form="formatted" )
   write(LDT_logunit,*) " -- Opening LDT Mask-Parameter Fill log file -- "

 end subroutine LDT_ParamMaskFill_Log


! =================
!BOP
!
! !ROUTINE: LDT_discreteParam_Fill
! \label{LDT_discreteParam_Fill}
!
! !INTERFACE:
  subroutine LDT_discreteParam_Fill ( n, pnc, pnr, gridtransform, pnl,  &
           param_value, undef, mask_value, fill_option, fill_value, fill_rad,&
           fill_value_antarctica, force_exclude_water)
! !USES:
    use ESMF
    use LDT_coreMod,   only : LDT_rc, LDT_config
    use LDT_fileIOMod, only : LDT_readDomainConfigSpecs
!
! !DESCRIPTION:
!
!  This routine ensures parameter-mask consistency with fill
!   options discrete-data parameter files.
!
! - Based on Yudong Tian's original code.
!
!EOP
    implicit none
!- Inputs:
    integer, intent(in)   :: n     ! Nest index
    integer, intent(in)   :: pnc   ! Parameter num of cols
    integer, intent(in)   :: pnr   ! Parameter num of rows
    integer, intent(in)   :: pnl   ! Parameter num of levels
    character(50), intent(in) :: gridtransform
    real,    intent(in)   :: undef
    real,    intent(in)   :: mask_value(pnc, pnr, pnl)
    real,    intent(inout):: param_value(pnc, pnr, pnl)
    real,    intent(in)   :: fill_value
    real,    intent(in)   :: fill_rad
    character(50), intent(in) :: fill_option
    real, optional, intent(in) :: fill_value_antarctica ! EMK
    logical, optional, intent(in) :: force_exclude_water ! EMK

!- Local parameters
    integer :: i, j, l
    integer :: i1, j1, ix, jx
    integer :: ifound, irad
    real    :: isum, isum2
    logical :: include_water ! EMK
! ______________________________________________________________

   write(fill_logunit,*) " -- Filling in discrete dataset's missing values --"
   write(fill_logunit,*) " -- Number of data layers: ", pnl
   write(fill_logunit,*) " -- Missing/water value: ", undef

!- Loop over output domain:
   do j = 1, pnr
      do i = 1, pnc

      !- Single layer discrete file:
         if( pnl == 1 ) then 

        != Mask: indicates land point;  Parameter: indicates "no surface type"
           if( mask_value(i, j, 1) == 1 .and. &
               param_value(i, j, 1) == undef )  then
              write(fill_logunit,*) "value needs filling: ", i, j

           !- Nearest neighbor fill option:
              if( fill_option == "neighbor" ) then
                 ifound = 0
                 do irad = 1, nint(fill_rad)
                    do j1 = j-irad, j+irad  ! search values N/S
                       jx = j1
                    !- Account for map edges:
                       if( jx <= 0   )  jx = 1
                       if( jx >= pnr )  jx = pnr
                       do i1 = i-irad, i+irad  ! search values E/W
                          ix = i1
                       !- Account for map edges:
                          if( ix <= 0  )  ix = pnc + ix
                          if( ix > pnc )  ix = ix - pnc
                          if( param_value(ix,jx,1) .ne. undef ) then 
                            ifound = ifound + 1
                            param_value(i, j, 1) = param_value(ix, jx, 1)
                            write(fill_logunit,*) "... found neighbor-value: ", param_value(i,j,1)
                            exit
                          end if
                       end do
                       if( ifound > 0 ) exit
                    end do
                 end do
                 if( ifound == 0 ) then
                    write(fill_logunit,*) " ... neighbor-value missing ... filling: "
                    param_value(i, j, 1) = fill_value
                end if
  
              else  ! Other fill options not currently supported
                write(LDT_logunit,*)" ERR:  For discrete-based parameter files, only the"
                write(LDT_logunit,*)"       nearest 'neighbor' option is currently available."
                write(LDT_logunit,*)" Stopping ..."
                call LDT_endrun
              end if   ! End fill option check

           end if

         !- N-layer file:  Tiled (discrete) data
        elseif( pnl > 1 ) then 

           ! EMK...Revise logic for excluding water points.  Allow exclusion
           ! for soil texture even if lake or water surface model type is
           ! used.
           include_water = LDT_rc%inc_water_pts
           if (present(force_exclude_water)) then
              if (force_exclude_water) then
                 include_water = .false.
              end if
           end if
                      
           isum = 0.
           if( include_water ) then
              !- INCLUDING WATER PIXELS
              isum = sum( param_value(i,j,1:pnl), &
                   mask=param_value(i,j,1:pnl).ne.LDT_rc%udef )
           else
              !- EXCLUDING WATER PIXELS
              do l = 1, pnl
                 if( l /= undef ) then
                    if( param_value(i,j,l) /= LDT_rc%udef ) then
                       isum = isum + param_value(i,j,l)
                    end if
                 end if
              end do
           endif
           
        != Mask: indicates land point;  Parameter: indicates "no surface type"
           if( mask_value(i, j, 1) == 1 .and. &
               isum < (1.0 - LDT_rc%gridcell_water_frac(n)) ) then
              
            select case ( fill_option )

           !- Nearest neighbor option:
              case( "neighbor" )
 
               do irad = 1, nint(fill_rad)
                  do j1 = j-irad, j+irad  ! search values N/S
                     jx = j1
                  !- Account for map edges:
                     if( jx <= 0 )   jx = 1
                     if( jx >= pnr ) jx = pnr
                     do i1 = i-irad, i+irad  ! search values E/W
                        ix = i1
                     !- Account for map edges:
                        if( ix <= 0 )  ix = pnc + ix
                        if( ix > pnc ) ix = ix - pnc
                     !- First neighboring pixel found with valid value array:
                        isum2 = 0
                        do l = 1, pnl
                           if( l /= undef ) then
                               if( param_value(ix,jx,l) /= LDT_rc%udef ) then
                                  isum2 = isum2 + param_value(ix,jx,l)
                               end if
                           end if
                        end do
                        if( undef > 0 .and. &
                            isum2>=(1.0 - LDT_rc%gridcell_water_frac(n)) ) then
!                            param_value(ix, jx, nint(undef)) < isum ) then
!                            param_value(ix, jx, nint(undef)) < 0.99 ) then
                            param_value(i, j, :) = param_value(ix, jx, :)
                            goto 199
                        elseif( undef <= 0 ) then
                            write(LDT_logunit,*) "(ParamMaskFill): 'neighbor' option "
                            write(LDT_logunit,*) " error ... reached a layer that is negative."
                            write(LDT_logunit,*) " Stopping .... "
                            call LDT_endrun
                            param_value(i, j, 1) = param_value(ix, jx, 1)
                            goto 199
                        end if
                     end do
                  end do
               end do

               Write(fill_logunit, *) "Can not find neighbor ! i= ", i, "  j=", j

             ! Since no neighbor found, fill with some user-specified value:
               param_value(i, j, :) = 0.0
               if( undef <= 0 ) then
                  param_value(i, j, 1) = 1.0
               else
                  ! EMK:  Allow optional different fill value south of 60S 
                  if (present(fill_value_antarctica)) then
                     if (LDT_domain(n)%lat(i + (j-1)*LDT_rc%lnc(n)) < -60.) then
                        param_value(i, j, nint(fill_value_antarctica)) = 1.0 
                     else
                        param_value(i, j, nint(fill_value)) = 1.0
                     end if
                  else
                     param_value(i, j, nint(fill_value)) = 1.0
                  end if
               endif

 199           continue
             
             case default   ! Other fill options not currently supported
                write(LDT_logunit,*)" ERR:  For discrete-based parameter files, only the"
                write(LDT_logunit,*)"       nearest 'neighbor' option is currently available."
                write(LDT_logunit,*)" Stopping ..."
                call LDT_endrun

             end select  ! End fill approach option

           != Mask == 0: Indicates water presence; 
           elseif( mask_value(i, j, 1) == 0 ) then
              param_value(i, j, :) = 0.0
              if( undef <= 0 ) then
                 param_value(i, j, 1) = 1.0
              else
                 param_value(i, j, nint(undef)) = 1.0
              endif
           end if

         end if  ! End num_layers condition
      end do     ! End output file loops
   end do

!== Mask: indicates water point;  Parameter: indicates valid value
   if( pnl == 1 ) then
     do j = 1, pnr
        do i = 1, pnc
          if( mask_value(i, j, 1) == 0 .and. &
             param_value(i, j, 1) /= undef )  then
             param_value(i, j, 1) = undef
!             param_value(i, j, 1) = LDT_rc%udef
          end if  ! End mask-param value check
        end do
     end do
   endif

 end subroutine LDT_discreteParam_Fill

! =================
!BOP
!
! !ROUTINE: LDT_contIndivParam_Fill
! \label{LDT_contIndivParam_Fill}
!
! !INTERFACE:
  subroutine LDT_contIndivParam_Fill( n, pnc, pnr, gridtransform, pnl,  &
       param_value, undef, mask_value, fill_option, fill_value, fill_rad, leave_good_data)

! !USES:
    use ESMF
    use LDT_coreMod,   only : LDT_rc, LDT_config
    use LDT_fileIOMod, only : LDT_readDomainConfigSpecs
!
! !DESCRIPTION:
!
!  This routine ensures parameter-mask consistency with fill
!   options for continuous data fields that are independent
!   of other parameters.
!
!EOP
    implicit none
!- Inputs:
    integer, intent(in)   :: n     ! Nest index
    integer, intent(in)   :: pnc   ! Parameter num of cols
    integer, intent(in)   :: pnr   ! Parameter num of rows
    integer, intent(in)   :: pnl   ! Parameter num of levels
    character(50), intent(in) :: gridtransform
    real,    intent(in)   :: undef
    real,    intent(in)   :: mask_value(pnc, pnr, pnl)
    real,    intent(inout):: param_value(pnc, pnr, pnl)
    real,    intent(in)   :: fill_value
    real,    intent(in)   :: fill_rad
    character(50), intent(in) :: fill_option
    logical, intent(in), optional :: leave_good_data

!- Local parameters:
    integer :: i, j, l
    integer :: i1, j1, ix, jx
    integer :: irad
    real    :: temp
    real    :: ifound(pnc, pnr)
    real    :: isum(pnc, pnr)
! ______________________________________________________________

   write(fill_logunit,*) " -- Filling in independent continuous dataset's missing values --"
   write(fill_logunit,*) " -- Number of data layers: ", pnl
   write(fill_logunit,*) " -- Missing/water value: ", undef

!- Loop over all available layers:
   do l = 1, pnl

      ifound = 0.
      isum   = 0.
  !== Mask: indicates land point;  Parameter: indicates "no surface type"
      do j = 1, pnr
         do i = 1, pnc

         != Skip here -- Mask: indicates water point;  Parameter: indicates valid value
            if( mask_value(i, j, 1) == 0 .and. &
                param_value(i, j, l) /= undef )  cycle

         != Mask: indicates land point;  Parameter: indicates "no surface type"
            if( mask_value(i, j, 1) == 1 .and. &
                param_value(i, j, l) == undef )  then
               write(fill_logunit,*) "value needs filling: ", l, i, j
 
            !- Average fill option:
               if( fill_option == "average" ) then
                  do irad = 1, nint(fill_rad)
                     do j1 = j-irad, j+irad   ! search values N/S
                        jx = j1
                     !- Account for map edges:
                        if( jx <= 0   ) jx = 1
                        if( jx >= pnr ) jx = pnr
                        do i1 = i-irad, i+irad  ! search values E/W
                           ix = i1
                        !- Account for map edges:
                           if( ix <= 0   ) ix = pnc + ix
                           if( ix > pnc ) ix = ix - pnc
                        !- Valid neighboring value found:
                           if( param_value(ix, jx, l) .ne. undef ) then
                              ifound(i,j) = ifound(i,j) + 1
                              isum(i,j) = isum(i,j) + param_value(ix, jx, l)
                           end if
                        end do
                     end do
                     if( ifound(i,j) >= 1 ) exit
                  end do
                  if( ifound(i,j) >= 1 ) then
                     param_value(i, j, l) = isum(i,j) / ifound(i,j)
                     write(fill_logunit,*)"filling value with average:",i,j,param_value(i, j, l)
                  end if

            !- Nearest neighbor option:
               elseif( fill_option == "neighbor" ) then
                  do irad = 1, nint(fill_rad)
                     do j1 = j-irad, j+irad   ! search values N/S
                        jx = j1
                     !- Account for map edges:
                        if( jx <= 0   ) jx = 1
                        if( jx >= pnr ) jx = pnr
                        do i1 = i-irad, i+irad  ! search values E/W
                           ix = i1
                        !- Account for map edges:
                           if( ix <= 0  ) ix = pnc + ix
                           if( ix > pnc ) ix = ix - pnc
                        !- Valid neighboring value found:
                           if( param_value(ix, jx, l) .NE. undef ) then
                              ifound(i,j) = ifound(i,j) + 1
                              param_value(i, j, l) = param_value(ix, jx, l)
                              write(fill_logunit,*)"filling value with neighbor:",&
                                    i,j,param_value(i,j,l)
                              exit
                           end if
                        end do
                     end do
                     if( ifound(i,j) >= 1 ) exit
                  end do  ! End radius search
               end if     ! End fill method option

            end if  ! End mask-param value check

         end do     ! End nc-loop
      end do        ! End nr-loop

  !== Ensure final agreement between mask and parameter values:
      do j = 1, pnr
         do i = 1, pnc

         !- Fill values where there are still inconsistencies:
            if( mask_value(i, j, 1) == 1 .and. &
                param_value(i, j, l) == undef )  then
               if( ifound(i,j) < 1. ) then
                  param_value(i, j, l) = fill_value
                  write(fill_logunit,*)"filling value with arb value:",&
                        i,j,param_value(i,j,l)
               end if
            end if  

            ! MMF scheme is active on lake/water grid cells. Thus, when MMF parameters are being processed,
            !  "leave_good_data" optional argument is used to leave acceptable parameter values in lake/water grid cells intact. 
            if(.not.present(leave_good_data)) then
               !== Mask: indicates water point;  Parameter: indicates valid value
               if( mask_value(i, j, 1) == 0 .and. &
                    param_value(i, j, l) /= undef )  then
                  param_value(i, j, l) = undef
               end if  ! End mask-param value check
            endif
         end do     ! End nc-loop
      end do        ! End nr-loop

   end do     ! End n-layers loop

 end subroutine LDT_contIndivParam_Fill

! ====================================
!
! !ROUTINE: LDT_contTileParam_Fill
! \label{LDT_contTileParam_Fill}
!
! !INTERFACE:
  subroutine LDT_contTileParam_Fill( n, pnc, pnr, &
           gridtransform, pnl, param_value, &
           param_fgrd, undef, mask_value, &
           fill_option, fill_value, fill_rad )
! !USES:
    use ESMF
    use LDT_coreMod,   only : LDT_rc, LDT_config
    use LDT_fileIOMod, only : LDT_readDomainConfigSpecs
!
! !DESCRIPTION:
!
!  This routine ensures parameter-mask consistency with fill
!   options for continuous data fields that are tiled and 
!   have co-dependent parameters (main tiled field+gridcell fractions).
!
!EOP
    implicit none
!- Inputs:
    integer, intent(in)   :: n     ! Nest index
    integer, intent(in)   :: pnc   ! Parameter num of cols
    integer, intent(in)   :: pnr   ! Parameter num of rows
    integer, intent(in)   :: pnl   ! Parameter num of levels
    character(50), intent(in) :: gridtransform
    real,    intent(in)   :: undef
    real,    intent(in)   :: mask_value(pnc, pnr, pnl)
    real,    intent(inout):: param_value(pnc, pnr, pnl)
    real,    intent(inout):: param_fgrd(pnc, pnr, pnl)
    real,    intent(in)   :: fill_value
    real,    intent(in)   :: fill_rad
    character(50), intent(in) :: fill_option

!- Local parameters:
    integer :: i, j, l
    integer :: i1, j1, ix, jx
    integer :: irad
    real    :: isum, isum_neighbor
! ______________________________________________________________

   write(fill_logunit,*) " -- Filling in tiled continuous dataset's missing values --"
   write(fill_logunit,*) " -- Number of data layers: ", pnl
   write(fill_logunit,*) " -- Missing/water value: ", undef

!- Loop over output points:
   do j = 1, pnr
      do i = 1, pnc

  ! ---  Tiled output data fields:
         if( gridtransform == "tile" ) then
            isum = 0
            isum = sum( param_fgrd(i,j,1:pnl), &
                      mask=param_fgrd(i,j,1:pnl).ne.undef )

         != Mask: indicates land point;  Parameter: indicates "no surface type"
            if( mask_value(i, j, 1) == 1 .and. &
                isum == 0 )  then
               write(fill_logunit,*) "value needs filling: ", i, j

           !- Nearest neighbor option:
              if( fill_option == "neighbor" ) then
outer:          do irad = 1, nint(fill_rad)
                  do j1 = j-irad, j+irad  ! search values N/S
                     jx = j1
                  !- Account for map edges:
                     if( jx <= 0   )  jx = 1
                     if( jx >= pnr )  jx = pnr
                     do i1 = i-irad, i+irad  ! search values E/W
                        ix = i1
                     !- Account for map edges:
                        if( ix <= 0  )  ix = pnc + ix
                        if( ix > pnc )  ix = ix - pnc
                     !- Sum nearest pixel layer-values:
                        isum_neighbor = 0
                        isum_neighbor = sum( param_fgrd(ix,jx,1:pnl), &
                                  mask=param_fgrd(ix,jx,1:pnl).ne.undef )

                     !- First neighboring pixel found with valid value array:
                        if( isum_neighbor > 0.98 ) then
                           param_fgrd(i, j, :)  = param_fgrd(ix, jx, :)
                           param_value(i, j, :) = param_value(ix, jx, :)
                           write(fill_logunit,*) " ... found neighbor-value: ", param_value(i, j, :)
                           exit outer
                        end if
                     end do
                  end do
               end do outer
               if( isum_neighbor <= 0.98 ) then
                  write(fill_logunit,*) " ... neighbor-value missing ... filling: "
                  param_fgrd(i, j, 1)  = 1.0
                  param_value(i, j, 1) = fill_value
               end if

              end if    ! End fill approach option

         != Mask: indicates water point;  Parameter: indicates "surface type"
            elseif( mask_value(i, j, 1) == 0 .and. &
                isum > 0.98 )  then
               param_fgrd(i, j, :)  = 0.
               param_value(i, j, :) = undef

            end if      ! End mask-parameter check

         else
            write(LDT_logunit,*) " (ParamMaskFill):  This continuous field is not tiled -- "
            write(LDT_logunit,*) " Stopping ..."
            call LDT_endrun
         end if      ! End grid-spatial transform condition
  
      end do     ! End output field loop
   end do

 end subroutine LDT_contTileParam_Fill

! ====================================

  
end module LDT_paramMaskCheckMod

