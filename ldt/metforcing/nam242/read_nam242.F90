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
!BOP
! !ROUTINE: read_nam242
! \label{read_nam242}
!
! !REVISION HISTORY:
!     Sep 2012: NOHRSC/NOAA: Initial specification
!  12 Mar 2013: James Geiger: Committed into LDT
!  10 Apr 2013: James Geiger: Rewrote to use GRIBAPI
!
! !INTERFACE:
subroutine read_nam242(n, findex, order, name00, name03, name06, &
                       F06flag, ferror, try)
! !USES:  
  use LDT_coreMod,        only : LDT_rc, LDT_domain
  use LDT_metforcingMod,  only : LDT_forc
  use LDT_logMod,         only : LDT_logunit, LDT_endrun, LDT_verify
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use nam242_forcingMod,  only : nam242_struc

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS:
  integer, intent(in)          :: n
  integer, intent(in)          :: findex
  integer, intent(in)          :: order
  character(len=*), intent(in) :: name00
  character(len=*), intent(in) :: name03
  character(len=*), intent(in) :: name06
  logical, intent(in)          :: F06flag
  integer, intent(out)         :: ferror 
  integer, intent(inout)       :: try
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  NAM forecast datasets, transforms into 9 LDT forcing 
!  parameters and interpolates to the LDT domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    index of the forcing source
!  \item[order]
!    flag indicating which data to be read (order=1, read for the previous 
!    3hr bookend, order=2, read for the next 3 hr bookend)
!  \item[name00]
!    name of the instantaneous forecast file
!  \item[name03]
!    name of the 3 hour NAM forecast file
!  \item[name06]
!    name of the 6 hour NAM forecast file
!  \item[F06flag]
!    flag to indicate if 6hr forecast data is required for this interval
!  \item[ferror]
!    flag to indicate success of the call (=1 indicates success)
!  \item[try]
!    index of the tries (in case of missing data)
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_nam242](\ref{interp_nam242}) \newline
!    spatially interpolates a NAM variable
!  \end{description}
!EOP
!==== Local Variables=======================
  
  character(len=LDT_CONST_PATH_LEN) :: fname
  integer :: lenfname
  integer :: lennamfname
  character(len=2) :: initcode
  character(len=2) :: fcstcode
  integer :: iv,c,r,t,k
  integer :: nforce
  integer :: ftn,iret
  integer :: gridDesc(200)
  integer :: nnam
  logical*1, allocatable :: lb(:)
  logical :: file_exists
  real, allocatable :: f(:)
  real, dimension(LDT_rc%lnc(n), LDT_rc%lnr(n)) :: varfield
  real    :: namdata_i(10,LDT_rc%ngrid(n))
  real    :: namdata_a(10,LDT_rc%ngrid(n))
  real    :: namdata_a_f06(10,LDT_rc%ngrid(n))
  real    :: value
  !integer :: nvars

  integer, dimension(nam242_struc(n)%nmif) :: idisc      ! discipline number of message
  integer, dimension(nam242_struc(n)%nmif) :: icateg     ! category number of message
  integer, dimension(nam242_struc(n)%nmif) :: iparam     ! parameter number of message
  integer, dimension(nam242_struc(n)%nmif) :: ileveltype ! level type
  integer, dimension(nam242_struc(n)%nmif) :: ilevel     ! level number
  integer, dimension(nam242_struc(n)%nmif) :: ipdtno     ! product definition template number
  integer, dimension(nam242_struc(n)%nmif) :: ifcst      ! forecast time
  integer, dimension(nam242_struc(n)%nmif) :: ifcst00    ! forecast time at fcst00 time
  integer, dimension(nam242_struc(n)%nmif) :: ifcst03    ! forecast time at fcst03 time
  integer, dimension(nam242_struc(n)%nmif) :: ifcst06    ! forecast time at fcst06 time
  integer, dimension(nam242_struc(n)%nmif) :: ifcst09    ! forecast time at fcst06 time
  integer, dimension(nam242_struc(n)%nmif) :: ifcst12    ! forecast time at fcst06 time

  integer :: igrib
  real    :: missingValue
  integer :: idisc_val
  integer :: icateg_val
  integer :: iparam_val
  integer :: ileveltype_val
  integer :: ilevel_val
  integer :: ipdtno_val
  integer :: ifcst_val
  logical :: var_found
  integer :: var_index
  logical, dimension(9) :: var_status

!=== End Variable Definition =======================

#if (defined USE_GRIBAPI)
  lennamfname = 26
  ifcst = -9999

  namdata_i = LDT_rc%udef
  namdata_a = LDT_rc%udef
  namdata_a_f06 = LDT_rc%udef
  varfield = 0.0
  nnam = (nam242_struc(n)%nc*nam242_struc(n)%nr)

!--------------------------------------------------------------------------
! Set the GRIB parameter specifiers
!--------------------------------------------------------------------------
  nforce = nam242_struc(n)%nmif

  ferror = 1

  ! order      : TMP, SPFH, DSWRF, DLWRF, UGRD, VGRD, PRES, APCP, ACPCP
  ! ave or inst: inst inst  inst   inst   inst  inst  inst  ave   ave
  idisc     =   (/000,000,  000,   000,   000,  000,  000,  000,  000/)
  icateg    =   (/000,001,  004,   005,   002,  002,  003,  001,  001/)
  iparam    =   (/000,000,  192,   192,   002,  003,  000,  008,  010/)
  ileveltype=   (/103,103,  001,   001,   103,  103,  001,  001,  001/)
  ilevel    =   (/002,002,  000,   000,   010,  010,  000,  000,  000/)
  ipdtno    =   (/000,000,  000,   000,   000,  000,  000,  008,  008/)
  ifcst00   =   (/000,000,  000,   000,   000,  000,  000,  000,  000/)
  ifcst03   =   (/003,003,  003,   003,   003,  003,  003,  000,  000/)
  ifcst06   =   (/006,006,  006,   006,   006,  006,  006,  003,  003/)
  ifcst09   =   (/009,009,  009,   009,   009,  009,  009,  006,  006/)
  ifcst12   =   (/012,012,  012,   012,   012,  012,  012,  009,  009/)

  ! GRIB-API mappings
  ! idisc      <--> discipline
  ! icateg     <--> parameterCategory
  ! iparam     <--> parameterNumber
  ! ileveltype <--> typeOfFirstFixedSurface
  ! ilevel     <--> level
  ! ipdtno     <--> productDefinitionTemplateNumber
  ! ifcst      <--> forecastTime

!--------------------------------------------------------------------------
! if there's a problem then ferror is set to zero
!--------------------------------------------------------------------------

  allocate(lb(nam242_struc(n)%nc*nam242_struc(n)%nr))
  allocate(f(nam242_struc(n)%nc*nam242_struc(n)%nr))

! Note that grib_count_in_file is not returning the correct number
! of messages in a grib file.  So instead of using nvars to control
! the loop that processes all the grib messages, check for the end of file.

!--------------------------------------------------------------------------
! read instantaneous fields
!--------------------------------------------------------------------------

  fname = trim(name00)
  !print*,'name00 = ',fname
  inquire (file=fname, exist=file_exists)

  if (file_exists) then
     call grib_open_file(ftn,trim(fname),'r',iret)
     call grib_multi_support_on
     if ( iret == 0 ) then
        lenfname = len(trim(fname))
        initcode = fname(lenfname-lennamfname-2:lenfname-lennamfname)
        fcstcode = fname(lenfname-lennamfname+6:lenfname-lennamfname+8)
        if (fcstcode == '00') then
           ifcst(:) = ifcst00(:)
        elseif (fcstcode == '03') then
           ifcst(:) = ifcst03(:)
        elseif (fcstcode == '06') then
           ifcst(:) = ifcst06(:)
        elseif (fcstcode == '09') then
           ifcst(:) = ifcst09(:)
        elseif (fcstcode == '12') then
           ifcst(:) = ifcst12(:)
        else
           write(LDT_logunit,*) 'ERR: read_nam242: Incorrect forecast time', &
                                fcstcode
        endif
     else
        write(LDT_logunit,*) 'ERR: read_nam242: Could not open ',trim(fname)
        ferror = 0
        deallocate(f)
        deallocate(lb)
        return
     endif

     !call grib_count_in_file(ftn,nvars,iret)
     !call LDT_verify(iret, 'ERR: read_nam242: Could not get number of vars')

     var_status = .false. 
     !do k = 1, nvars
     do

        call grib_new_from_file(ftn, igrib, iret)
        if ( iret == GRIB_END_OF_FILE .or. all(var_status) ) then
           call grib_release(igrib,iret)
           exit
        endif

        if ( iret /= 0 ) then 
           write(LDT_logunit,*) &
                'ERR: read_nam242: Could not retrieve entries in file: ', &
                trim(fname)
           ferror = 0
           deallocate(lb)
           deallocate(f)
           return           
        endif

        call grib_get(igrib,'discipline',idisc_val,iret)
        call LDT_verify(iret,'ERR: read_nam242: could not get discipline value')

        call grib_get(igrib,'parameterCategory',icateg_val,iret)
        call LDT_verify(iret, &
                      'ERR: read_nam242: could not get parameterCategory value')

        call grib_get(igrib,'parameterNumber',iparam_val,iret)
        call LDT_verify(iret, &
                        'ERR: read_nam242: could not get parameterNumber value')

        call grib_get(igrib,'typeOfFirstFixedSurface',ileveltype_val,iret)
        call LDT_verify(iret, &
                'ERR: read_nam242: could not get typeOfFirstFixedSurface value')

        call grib_get(igrib,'level',ilevel_val,iret)
        call LDT_verify(iret, 'ERR: read_nam242: could not get level value')

        call grib_get(igrib,'productDefinitionTemplateNumber',ipdtno_val,iret)
        call LDT_verify(iret, &
        'ERR: read_nam242: could not get productDefinitionTemplateNumber value')

        call grib_get(igrib,'forecastTime',ifcst_val,iret)
        call LDT_verify(iret, 'ERR: read_nam242: could not get forecastTime value')


        var_found = .false. 
        do iv = 1, 9
           if ( ( idisc_val      == idisc(iv)      ) .and. &
                ( icateg_val     == icateg(iv)     ) .and. &
                ( iparam_val     == iparam(iv)     ) .and. &
                ( ileveltype_val == ileveltype(iv) ) .and. &
                ( ilevel_val     == ilevel(iv)     ) .and. &
                ( ipdtno_val     == ipdtno(iv)     ) .and. &
                ( ifcst_val      == ifcst(iv) )    ) then
              var_found = .true.
              var_index = iv
              var_status(iv) = .true. 
              exit
           endif
        enddo

        if ( var_found ) then
           f = -9999.0
           call grib_get(igrib,'values',f,iret)

           if(iret /= 0) then 
              write(LDT_logunit,*) &
                  'ERR: read_nam242: Could not retrieve values in file: ', &
                  trim(fname)
              write(LDT_logunit,*) &
                  'ERR: read_nam242: Could not retrieve values for forcing: ',&
                  var_index
              ferror = 0
              deallocate(lb)
              deallocate(f)
              return           
           endif

           call grib_get(igrib,'missingValue',missingValue,iret)
           if ( iret /= 0 ) then
           !write(LDT_logunit,*) 'ERR: read_nam242: Could not get missingValue'
              missingValue = LDT_rc%udef
           endif

           lb = .false.
           do t = 1, nnam
              if ( f(t) /= missingValue ) then
                 lb(t) = .true.
              endif
           enddo
!--------------------------------------------------------------------------
! If field successfully retrieved, interplate to LDT domain
!--------------------------------------------------------------------------
           call interp_nam242(n,findex,iparam(var_index),nnam,f,lb,&
                              LDT_rc%gridDesc(n,:),&
                              LDT_rc%lnc(n),LDT_rc%lnr(n),varfield)

           do r = 1, LDT_rc%lnr(n)
              do c = 1, LDT_rc%lnc(n)
                 if ( LDT_domain(n)%gindex(c,r) /= -1 ) then 
                    namdata_i(iv,LDT_domain(n)%gindex(c,r)) = varfield(c,r)
                 endif
              enddo
           enddo
        endif

        call grib_release(igrib,iret)
        call LDT_verify(iret, 'ERR: read_nam242: Could not release igrib')
     enddo

     call grib_close_file(ftn)

     do k = 1, 9
        if ( .not. var_status(k) ) then 
           write(LDT_logunit,*) &
             'ERR: read_nam242: Could not retrieve all inst entries in file: ',&
             trim(fname)
           write(LDT_logunit,*) &
                'ERR: read_nam242: retrieval status ',var_status
           ferror = 0
           deallocate(lb)
           deallocate(f)
           return
        endif
     enddo
  else
     write(LDT_logunit,*) 'ERR: read_nam242: ',trim(fname),' does not exist.'
     ferror = 0
     deallocate(f)
     deallocate(lb)
     return
  endif

!--------------------------------------------------------------------------
! read 3hr forecast for time averaged fields
!--------------------------------------------------------------------------

  fname = trim(name03)
  !print*,'name03 = ',fname
  inquire (file=fname, exist=file_exists)

  if (file_exists) then
     call grib_open_file(ftn,trim(fname),'r',iret)
     call grib_multi_support_on
     if ( iret == 0 ) then
        lenfname = len(trim(fname))
        initcode = fname(lenfname-lennamfname-2:lenfname-lennamfname)
        fcstcode = fname(lenfname-lennamfname+6:lenfname-lennamfname+8)
        if (fcstcode == '00') then
           ifcst(:) = ifcst00(:)
        elseif (fcstcode == '03') then
           ifcst(:) = ifcst03(:)
        elseif (fcstcode == '06') then
           ifcst(:) = ifcst06(:)
        elseif (fcstcode == '09') then
           ifcst(:) = ifcst09(:)
        elseif (fcstcode == '12') then
           ifcst(:) = ifcst12(:)
        else
           write(LDT_logunit,*) 'ERR: read_nam242: Incorrect forecast time', &
                                fcstcode
        endif
     else
        write(LDT_logunit,*) 'ERR: read_nam242: Could not open ',trim(fname)
        ferror = 0
        deallocate(f)
        deallocate(lb)
        return
     endif

     !call grib_count_in_file(ftn,nvars,iret)
     !call LDT_verify(iret, 'ERR: read_nam242: Could not get number of vars')

     var_status = .false. 
     !do k = 1, nvars
     do

        call grib_new_from_file(ftn, igrib, iret)
        if ( iret == GRIB_END_OF_FILE .or. all(var_status) ) then
           call grib_release(igrib,iret)
           exit
        endif

        if ( iret /= 0 ) then 
           write(LDT_logunit,*) &
           'ERR: read_nam242: Could not retrieve entries in file: ',trim(fname)
           ferror = 0
           deallocate(lb)
           deallocate(f)
           return           
        endif

        call grib_get(igrib,'discipline',idisc_val,iret)
        call LDT_verify(iret,'ERR: read_nam242: could not get discipline value')

        call grib_get(igrib,'parameterCategory',icateg_val,iret)
        call LDT_verify(iret, &
                      'ERR: read_nam242: could not get parameterCategory value')

        call grib_get(igrib,'parameterNumber',iparam_val,iret)
        call LDT_verify(iret, &
                        'ERR: read_nam242: could not get parameterNumber value')

        call grib_get(igrib,'typeOfFirstFixedSurface',ileveltype_val,iret)
        call LDT_verify(iret, &
                'ERR: read_nam242: could not get typeOfFirstFixedSurface value')

        call grib_get(igrib,'level',ilevel_val,iret)
        call LDT_verify(iret, 'ERR: read_nam242: could not get level value')

        call grib_get(igrib,'productDefinitionTemplateNumber',ipdtno_val,iret)
        call LDT_verify(iret, &
        'ERR: read_nam242: could not get productDefinitionTemplateNumber value')

        call grib_get(igrib,'forecastTime',ifcst_val,iret)
        call LDT_verify(iret, 'ERR: read_nam242: could not get forecastTime value')


        var_found = .false. 
        do iv = 1, 9
           if ( ( idisc_val      == idisc(iv)      ) .and. &
                ( icateg_val     == icateg(iv)     ) .and. &
                ( iparam_val     == iparam(iv)     ) .and. &
                ( ileveltype_val == ileveltype(iv) ) .and. &
                ( ilevel_val     == ilevel(iv)     ) .and. &
                ( ipdtno_val     == ipdtno(iv)     ) .and. &
                ( ifcst_val      == ifcst(iv) )    ) then
              var_found = .true.
              var_index = iv
              var_status(iv) = .true. 
              exit
           endif
        enddo

        if ( var_found ) then
           f = -9999.0
           call grib_get(igrib,'values',f,iret)

           if(iret /= 0) then 
              write(LDT_logunit,*) &
                  'ERR: read_nam242: Could not retrieve values in file: ',&
                  trim(fname)
              write(LDT_logunit,*) &
                  'ERR: read_nam242: Could not retrieve values for forcing: ',&
                  var_index
              ferror = 0
              deallocate(lb)
              deallocate(f)
              return           
           endif

           call grib_get(igrib,'missingValue',missingValue,iret)
           if ( iret /= 0 ) then
           !write(LDT_logunit,*) 'ERR: read_nam242: Could not get missingValue'
              missingValue = LDT_rc%udef
           endif

           lb = .false.
           do t = 1, nnam
              if ( f(t) /= missingValue ) then
                 lb(t) = .true.
              endif
           enddo
!--------------------------------------------------------------------------
! If field successfully retrieved, interplate to LDT domain
!--------------------------------------------------------------------------
           call interp_nam242(n,findex,iparam(var_index),nnam,f,lb,&
                              LDT_rc%gridDesc(n,:),&
                              LDT_rc%lnc(n),LDT_rc%lnr(n),varfield)

           do r = 1, LDT_rc%lnr(n)
              do c = 1, LDT_rc%lnc(n)
                 if ( LDT_domain(n)%gindex(c,r) /= -1 ) then 
                    namdata_a(iv,LDT_domain(n)%gindex(c,r)) = varfield(c,r)
                 endif
              enddo
           enddo
        endif

        call grib_release(igrib,iret)
        call LDT_verify(iret, 'ERR: read_nam242: Could not release igrib')
     enddo

     call grib_close_file(ftn)

     do k = 1, 9
        if ( .not. var_status(k) ) then 
           write(LDT_logunit,*) &
              'ERR: read_nam242: Could not retrieve all 3hr entries in file: ',&
              trim(fname)
           write(LDT_logunit,*) &
                'ERR: read_nam242: retrieval status ',var_status
           ferror = 0
           deallocate(lb)
           deallocate(f)
           return
        endif
     enddo
  else
     write(LDT_logunit,*) 'ERR: read_nam242: ',trim(fname),' does not exist.'
     ferror = 0
     deallocate(f)
     deallocate(lb)
     return
  endif

!--------------------------------------------------------------------------
! read 6hr forecast for time averaged fields, if required. 
!--------------------------------------------------------------------------

  if ( F06flag ) then 
   fname = trim(name06)
   !print*,'name06 = ',fname
   inquire (file=fname, exist=file_exists)

   if (file_exists) then
      call grib_open_file(ftn,trim(fname),'r',iret)
      call grib_multi_support_on
      if ( iret == 0 ) then
         lenfname = len(trim(fname))
         initcode = fname(lenfname-lennamfname-2:lenfname-lennamfname)
         fcstcode = fname(lenfname-lennamfname+6:lenfname-lennamfname+8)
         if (fcstcode == '00') then
            ifcst(:) = ifcst00(:)
         elseif (fcstcode == '03') then
            ifcst(:) = ifcst03(:)
         elseif (fcstcode == '06') then
            ifcst(:) = ifcst06(:)
         elseif (fcstcode == '09') then
            ifcst(:) = ifcst09(:)
         elseif (fcstcode == '12') then
            ifcst(:) = ifcst12(:)
         else
            write(LDT_logunit,*) 'ERR: read_nam242: Incorrect forecast time', &
                                 fcstcode
         endif
      else
         write(LDT_logunit,*) 'ERR: read_nam242: Could not open ',trim(fname)
         ferror = 0
         deallocate(f)
         deallocate(lb)
         return
      endif

      !call grib_count_in_file(ftn,nvars,iret)
      !call LDT_verify(iret, 'ERR: read_nam242: Could not get number of vars')

      var_status = .false. 
      !do k = 1, nvars
      do

         call grib_new_from_file(ftn, igrib, iret)
         if ( iret == GRIB_END_OF_FILE .or. all(var_status) ) then
            call grib_release(igrib,iret)
            exit
         endif

         if ( iret /= 0 ) then 
            write(LDT_logunit,*) &
            'ERR: read_nam242: Could not retrieve entries in file: ',trim(fname)
            ferror = 0
            deallocate(lb)
            deallocate(f)
            return           
         endif

         call grib_get(igrib,'discipline',idisc_val,iret)
         call LDT_verify(iret, 'ERR: read_nam242: could not get discipline value')

         call grib_get(igrib,'parameterCategory',icateg_val,iret)
         call LDT_verify(iret, &
                      'ERR: read_nam242: could not get parameterCategory value')

         call grib_get(igrib,'parameterNumber',iparam_val,iret)
         call LDT_verify(iret, &
                        'ERR: read_nam242: could not get parameterNumber value')

         call grib_get(igrib,'typeOfFirstFixedSurface',ileveltype_val,iret)
         call LDT_verify(iret, &
                'ERR: read_nam242: could not get typeOfFirstFixedSurface value')

         call grib_get(igrib,'level',ilevel_val,iret)
         call LDT_verify(iret, 'ERR: read_nam242: could not get level value')

         call grib_get(igrib,'productDefinitionTemplateNumber',ipdtno_val,iret)
         call LDT_verify(iret, &
         'ERR: read_nam242: could not get productDefinitionTemplateNumber value')

         call grib_get(igrib,'forecastTime',ifcst_val,iret)
         call LDT_verify(iret, 'ERR: read_nam242: could not get forecastTime value')


         var_found = .false. 
         do iv = 1, 9
            if ( ( idisc_val      == idisc(iv)      ) .and. &
                 ( icateg_val     == icateg(iv)     ) .and. &
                 ( iparam_val     == iparam(iv)     ) .and. &
                 ( ileveltype_val == ileveltype(iv) ) .and. &
                 ( ilevel_val     == ilevel(iv)     ) .and. &
                 ( ipdtno_val     == ipdtno(iv)     ) .and. &
                 ( ifcst_val      == ifcst(iv) )    ) then
               var_found = .true.
               var_index = iv
               var_status(iv) = .true. 
               exit
            endif
         enddo

         if ( var_found ) then
            f = -9999.0
            call grib_get(igrib,'values',f,iret)

            if(iret /= 0) then 
              write(LDT_logunit,*) &
                  'ERR: read_nam242: Could not retrieve values in file: ',&
                  trim(fname)
              write(LDT_logunit,*) &
                  'ERR: read_nam242: Could not retrieve values for forcing: ',&
                  var_index
               ferror = 0
               deallocate(lb)
               deallocate(f)
               return           
            endif

            call grib_get(igrib,'missingValue',missingValue,iret)
            if ( iret /= 0 ) then
            !write(LDT_logunit,*) 'ERR: read_nam242: Could not get missingValue'
               missingValue = LDT_rc%udef
            endif

            lb = .false.
            do t = 1, nnam
               if ( f(t) /= missingValue ) then
                  lb(t) = .true.
               endif
            enddo
 !--------------------------------------------------------------------------
 ! If field successfully retrieved, interplate to LDT domain
 !--------------------------------------------------------------------------
            call interp_nam242(n,findex,iparam(var_index),nnam,f,lb,&
                               LDT_rc%gridDesc(n,:),&
                               LDT_rc%lnc(n),LDT_rc%lnr(n),varfield)

            do r = 1, LDT_rc%lnr(n)
               do c = 1, LDT_rc%lnc(n)
                  if ( LDT_domain(n)%gindex(c,r) /= -1 ) then 
                     namdata_a_f06(iv,LDT_domain(n)%gindex(c,r)) = varfield(c,r)
                  endif
               enddo
            enddo
         endif

         call grib_release(igrib,iret)
         call LDT_verify(iret, 'ERR: read_nam242: Could not release igrib')
      enddo

      call grib_close_file(ftn)

      do k = 1, 9
         if ( .not. var_status(k) ) then 
            write(LDT_logunit,*) &
              'ERR: read_nam242: Could not retrieve all 6hr entries in file: ',&
              trim(fname)
            write(LDT_logunit,*) &
                 'ERR: read_nam242: retrieval status ',var_status
            ferror = 0
            deallocate(lb)
            deallocate(f)
            return
         endif
      enddo
   else
      write(LDT_logunit,*) 'ERR: read_nam242: ',trim(fname),' does not exist.'
      ferror = 0
      deallocate(f)
      deallocate(lb)
      return
   endif
  endif

  deallocate(f)
  deallocate(lb)

!--------------------------------------------------------------------------
! Place the interpolated data into the LDT arrays
!--------------------------------------------------------------------------

  do iv=1,nam242_struc(n)%nmif
     do t=1,LDT_rc%ngrid(n)
        if ( F06flag ) then 
           if(iv.eq.8.or.iv.eq.9) then ! 3-6 hour time avgd: APCP, ACPCP
              if ( namdata_a_F06(iv,t) /= LDT_rc%udef ) then
                 value = namdata_a_F06(iv,t)/10800.0
              else
                 value = LDT_rc%udef
              endif
              if ( order == 1 ) then
                 LDT_forc(n,findex)%metdata1(iv,t) = value
              else
                 LDT_forc(n,findex)%metdata2(iv,t) = value
              endif
           else ! instantaneous
              if ( order == 1 ) then
                 LDT_forc(n,findex)%metdata1(iv,t) = namdata_i(iv,t)
              else
                 LDT_forc(n,findex)%metdata2(iv,t) = namdata_i(iv,t)
              endif
           endif
        else
           if(iv.eq.8.or.iv.eq.9) then ! 0-3 hour time avgd: APCP, ACPCP
              if ( namdata_a(iv,t) /= LDT_rc%udef ) then
                 value = namdata_a(iv,t)/10800.0
              else
                 value = LDT_rc%udef
              endif
              if ( order == 1 ) then
                 LDT_forc(n,findex)%metdata1(iv,t) = value
              else
                 LDT_forc(n,findex)%metdata2(iv,t) = value
              endif
           else
              if ( order == 1 ) then
                 LDT_forc(n,findex)%metdata1(iv,t) = namdata_i(iv,t)
              else
                 LDT_forc(n,findex)%metdata2(iv,t) = namdata_i(iv,t)
              endif
           endif
        endif
     enddo
  enddo

#else
  write(LDT_logunit,*) 'ERR: read_nam242 requires GRIB-API'
  write(LDT_logunit,*) 'ERR: please recompile LDT'
  call LDT_endrun
#endif

end subroutine read_nam242


!BOP
! !ROUTINE: interp_nam242
! \label{interp_nam242}
!
! !INTERFACE:
subroutine interp_nam242(n,findex,iparam,nnam,f,lb,ldt_gds,nc,nr,varfield)
! !USES:
  use LDT_coreMod,       only : LDT_rc, LDT_domain
  use nam242_forcingMod, only : nam242_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
  integer, intent(in) :: iparam
  integer, intent(in) :: nnam
  real, intent(out)   :: f(nnam)
  logical*1           :: lb(nnam)
  real                :: ldt_gds(20)
  integer, intent(in) :: nc
  integer, intent(in) :: nr
  real, intent(out)   :: varfield(nc,nr)
!
! !DESCRIPTION:
!   This subroutine interpolates a given NAM field 
!   to the LDT grid. 
!  The arguments are: 
!  \begin{description}
! \item[n]
!  index of the nest
! \item[n]
!  index of the forcing source
! \item[iparam]
!  paramter indicator
! \item[nnam]
!  number of elements in the input grid
! \item[f]
!  input data array to be interpolated
! \item[lb]
!  input bitmap
! \item[ldt\_gds]
!  array description of the LDT grid
! \item[nc]
!  number of columns (in the east-west dimension) in the LDT grid
! \item[nr]
!  number of rows (in the north-south dimension) in the LDT grid
! \item[varfield]
!  output interpolated field
!  \end{description} 
! 
!
!  The routines invoked are: 
!  \begin{description}
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
! \end{description}
!EOP
  integer :: iret
  integer :: count1,i,j,mo

  real, dimension(nc*nr) :: ldt1d

  logical*1 :: lo(nc*nr)

!=== End variable declarations
  mo = nc*nr

!-----------------------------------------------------------------------
! Initialize output bitmap. Important for soil moisture and temp.
!-----------------------------------------------------------------------
  lo = .true.

!-----------------------------------------------------------------------  
! Interpolate to LDT grid
!-----------------------------------------------------------------------  
  if(LDT_rc%met_gridtransform(findex).eq."bilinear") then 
     call bilinear_interp(ldt_gds,lb,f,lo,ldt1d,nam242_struc(n)%mi,mo,&
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          nam242_struc(n)%w111,nam242_struc(n)%w121,&
          nam242_struc(n)%w211,nam242_struc(n)%w221,&
          nam242_struc(n)%n111,nam242_struc(n)%n121,&
          nam242_struc(n)%n211,nam242_struc(n)%n221,LDT_rc%udef, iret)
  elseif(LDT_rc%met_gridtransform(findex).eq."budget-bilinear") then 
     if (iparam==7 .or. iparam==196)then     
        call conserv_interp(ldt_gds,lb,f,lo,ldt1d,nam242_struc(n)%mi,mo, & 
             LDT_domain(n)%lat, LDT_domain(n)%lon,&
             nam242_struc(n)%w112,nam242_struc(n)%w122,&
             nam242_struc(n)%w212,nam242_struc(n)%w222,&
             nam242_struc(n)%n112,nam242_struc(n)%n122,&
             nam242_struc(n)%n212,nam242_struc(n)%n222,LDT_rc%udef,iret)
     else 
        call bilinear_interp(ldt_gds,lb,f,lo,ldt1d,nam242_struc(n)%mi,mo,&
             LDT_domain(n)%lat, LDT_domain(n)%lon,&
             nam242_struc(n)%w111,nam242_struc(n)%w121,&
             nam242_struc(n)%w211,nam242_struc(n)%w221,&
             nam242_struc(n)%n111,nam242_struc(n)%n121,&
             nam242_struc(n)%n211,nam242_struc(n)%n221,LDT_rc%udef,iret)
     endif
  endif

!-----------------------------------------------------------------------    
! Create 2D array for main program. Also define a "soil" mask
! due to different geography between NAM & LDAS. For LDAS land 
! points not included in NAM geography dataset only.
!-----------------------------------------------------------------------    
  count1 = 0
  do j = 1, nr
     do i = 1, nc
        varfield(i,j) = ldt1d(i+count1)
     enddo
     count1 = count1 + nc
  enddo

end subroutine interp_nam242
