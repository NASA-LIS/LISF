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
!
! !ROUTINE: read_ecmwf 
!  \label{read_ecmwf}
!
! !REVISION HISTORY:
!  18 Jun 2003: Urszula Jambor; original code
!  22 Jun 2011: Hiroko Beaudoing; Updated to GRIB_API routines 
!               to handle both GRIB1 and GRIB2 data
!
! !INTERFACE:
subroutine read_ecmwf( order, n, findex, &
     avgfile1, avgfile2, instfile, yr, mon, da, hr, ferror )
! !USES:
  use LDT_coreMod, only        : LDT_rc, LDT_domain
  use LDT_logMod, only         : LDT_logunit, LDT_endrun
  use LDT_metforcingMod, only : LDT_forc
  use ecmwf_forcingMod,   only : ecmwf_struc

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: order     ! lower(1) or upper(2) time interval bdry
  integer, intent(in)    :: n         ! lower(1) or upper(2) time interval bdry
  integer, intent(in)    :: findex
  integer, intent(in)    :: yr,mon,da,hr ! data and hour (multiple of 3)
  integer, intent(inout) :: ferror    ! set to zero if there's an error
  character(len=*)       :: avgfile1
  character(len=*)       :: avgfile2
  character(len=*)       :: instfile

! !DESCRIPTION:
!  For the given time, reads the parameters from 1/4 degree
!  ECMWF data, transforms into 9 LDT forcing 
!  parameters and interpolates to the LDT domain.

!  ECMWF model output variables used to force LDT fall into 2 categories:
!  1) inst3, Instantaneous values, available every 3 hours \newline
!  2) accum, Time integrated values (accumulations), updated every 3 hours \newline
!
!  NOTE 1: be aware that ECMWF outputs large-scale and convective precipitation
!  separately.  For total precipitation, need to sum the two fields,
!  LSP+CP=TP. \newline
!  NOTE 2: only time2 SW flux accumulations used in interpolation \newline
!  NOTE 3: read in of Albedo is currently suppressed.  This field is 
!  instantaneous and available every 6 hours.  At this time, all LDT LSMs 
!  replace this parameter. \newline
!
!  ECMWF FORCING VARIABLES: \newline
!  1. T         inst3, near-surfcae temperature, ~10 metres [K] \newline
!  2. Q         inst3, near-surface specific humidity, ~10 metres[kg/kg] \newline
!  3. SSRD      accum, surface solar radiation downwards [W m**-2 s] \newline
!  4. STRD      accum, surface thermal radiation downwards [W m**-2 s] \newline
!  5. U         inst3, zonal wind,~10 metres [m/s] \newline
!  6. V         inst3, meridional wind,~10 metres[m/s] \newline
!  7. SP        inst3, surface pressure [Pa]  \newline
!  8. LSP       accum, large scale precipitation [m] \newline
!  9. CP        accum, convective precipitation [m] \newline
!
!
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous 
!    3 hourly instance, order=2, read the next 3 hourly instance)
!  \item[n]
!    index of the nest
!  \item[yr]
!    current year
!  \item[mon]
!    current month
!  \item[da]
!    current day of the year
!  \item[hr]
!    current hour of day
!  \item[ferror]
!    flag to indicate success of the call (=0 indicates success)
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[retrieve\_inst\_ecmwfvars](\ref{retrieve_inst_ecmwfvars}) \newline
!    retrieves the instantaneous variable
!  \item[retrieve\_accum\_ecmwfvars](\ref{retrieve_accum_ecmwfvars}) \newline
!    retrieves the accumulated variable
!  \item[interp\_ecmwf](\ref{interp_ecmwf}) \newline
!    spatial interpolation from the ECMWF grid to the LDT grid 
!  \end{description}
!EOP
  integer                         :: c,r
  integer                         :: iv,t
  integer                         :: necmwf
  integer                         :: fret,ferror1,ferror2
  real                            :: glbdata_i(9,LDT_rc%ngrid(n))
  real                            :: glbdata_a(9,LDT_rc%ngrid(n))

  necmwf = ecmwf_struc(n)%ncold*ecmwf_struc(n)%nrold
  
  glbdata_i = LDT_rc%udef
  glbdata_a = LDT_rc%udef

  ferror1 = 1
  ferror2 = 1
  fret = 0

  call retrieve_inst_ecmwfvars(n, findex, instfile, glbdata_i, fret)
  if (fret.ne.0) then
     ferror1 = 0
  endif

  if ( order.eq.2 ) then
  call retrieve_accum_ecmwfvars(n, findex, avgfile1, avgfile2, glbdata_a, &
                                hr, fret)
   if (fret.ne.0) then
      ferror2 = 0
   endif
  endif  ! order==2

  if ( ferror1 == 1 .and. ferror2 == 1) then
   ferror = 1     ! success
  else
   ferror = 0     ! problem, roll back one day
  endif
  if ( ferror == 1 ) then  ! only proceed if retrieve calls were successful
   do iv=1,9
     do t=1,LDT_rc%ngrid(n)
        ! these are time avgd fields
        if ( iv.eq.3 .or. iv.eq.4 .or. iv.eq.8 .or. iv.eq.9 ) then
           if(order.eq.1) then 
              LDT_forc(n,findex)%metdata1(iv,t) = glbdata_a(iv,t)
           else
              LDT_forc(n,findex)%metdata2(iv,t) = glbdata_a(iv,t)
           endif
        ! these are instantaneous
        else
           if(order.eq.1) then 
              LDT_forc(n,findex)%metdata1(iv,t) = glbdata_i(iv,t)
           else
              LDT_forc(n,findex)%metdata2(iv,t) = glbdata_i(iv,t)
           endif
        endif
     enddo
   enddo
  endif    ! ferror == 1

end subroutine read_ecmwf
  
!BOP
! !ROUTINE: retrieve_inst_ecmwfvars
! \label{retrieve_inst_ecmwfvars}
!
! !INTERFACE:
subroutine retrieve_inst_ecmwfvars(n, findex, instfile, glbdata, fret)
  use LDT_coreMod,      only : LDT_rc, LDT_domain
  use LDT_logMod,       only : LDT_logunit, LDT_verify
  use ecmwf_forcingMod, only : ecmwf_struc
#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS:
  integer,   intent(in)        :: n 
  integer,   intent(in)        :: findex
  character(len=*), intent(in) :: instfile
  real                         :: glbdata(9,LDT_rc%ngrid(n))
  integer, intent(inout)       :: fret
!
! !DESCRIPTION: 
! This routine opens the corresponding ECMWF data file to extract
! the specified variable, which represents an instantaneous value.
! Should be used for near-surface temperature, specific humidity,
! winds, and surface pressure.
! 
!EOP
  integer,parameter  :: N_INST_VARS=5
  integer            :: necmwf
  integer            :: iret
  integer            :: c,r,iv,ftn,igrib
  integer            :: kk,nvars
  character(len=20)  :: shortName 
  real               :: missingValue 
  logical            :: file_exists
  character(len=2)   :: varname(N_INST_VARS)
  logical            :: var_status(N_INST_VARS)
  real,  allocatable     :: f(:)
  logical*1, allocatable :: lb(:)
  integer            :: var_index
  logical            :: var_found
  logical            :: pcp_flag
  integer            :: rel_index(N_INST_VARS)
  real, dimension(LDT_rc%lnc(n), LDT_rc%lnr(n))  :: varfield

  ! Then instantaneous files contain only 6 GRIB messages
  ! 1 = t
  ! 2 = q
  ! 3 = u
  ! 4 = v
  ! 5 = sp
  ! 6 = al <-- not used by LDT

  !=== set GRIB shortName 
  varname = (/ "t ","q ","u ","v ","sp" /)
  rel_index = (/ 1, 2, 5, 6, 7 /) ! index of variable w.r.t. the
                                  ! full list of 9 forcing variables

  necmwf = ecmwf_struc(n)%ncold*ecmwf_struc(n)%nrold

  varfield = 0
  var_status = .false.

#if (defined USE_GRIBAPI)
  !=== Set up to open file and retrieve specified field 
  fret = -1
  inquire(file=instfile,exist=file_exists) 
  if ( file_exists ) then 
     call grib_open_file(ftn,instfile,'r',iret)
     call LDT_verify(iret, 'failed to open file '//trim(instfile))

     call grib_count_in_file(ftn,nvars,iret)
     call LDT_verify(iret, 'error in grib_count_in_file in retrieve_inst_ecmwfvars')
     
     allocate(lb(necmwf))
     allocate(f(necmwf))

     do kk=1,nvars
        call grib_new_from_file(ftn,igrib,iret)
        call grib_get(igrib,'shortName',shortName)
        
        var_found = .false. 
        do iv=1,N_INST_VARS            
           if ( shortName == varname(iv) ) then
              var_found = .true. 
              var_index = rel_index(iv)
              var_status(iv) = .true. 
              exit
           endif
        enddo
        
        f = -9999.0
        call grib_get(igrib,'values',f,iret)
        call LDT_verify(iret, 'failed to get values in retrieve_inst_ecmwfvars')

        if ( iret .ne. 0 ) then   ! read not successful
         write(LDT_logunit,*) 'problem reading inst for ',shortName,iret
         call grib_release(igrib)
         deallocate(lb)
         deallocate(f)     
         fret = -1
         return 
        endif
 
        call grib_get(igrib,'missingValue',missingValue,iret)
        call LDT_verify(iret,'failed to get missingValue in retrieve_inst_ecmwfvars')
        call grib_release(igrib)
        call LDT_verify(iret,'error in grib_release in read_ecmwf')

        if(var_found) then
           lb = .false. 
           do c = 1, necmwf
              if ( f(c) .ne. missingValue ) then
                 lb(c) = .true.
              endif
           end do
           
           pcp_flag = .false.
           
           call interp_ecmwf(n,findex,pcp_flag,necmwf,f, &
                             lb,LDT_rc%gridDesc(n,:),    &
                             LDT_rc%lnc(n),LDT_rc%lnr(n),varfield)
           
           do r=1,LDT_rc%lnr(n)
              do c=1,LDT_rc%lnc(n)
                 if ( LDT_domain(n)%gindex(c,r) /= -1 ) then 
                    glbdata(var_index,LDT_domain(n)%gindex(c,r)) = &
                       varfield(c,r)
                 endif
              enddo
           enddo
        endif   
     enddo
     call grib_close_file(ftn)
  
     deallocate(lb)
     deallocate(f)     
     
     do kk=1,N_INST_VARS
        if ( .not. var_status(kk) ) then 
           write(LDT_logunit,*) &
                'ERR: Could not retrieve all entries in file: ',trim(instfile)
           fret = -1
           return
        endif
     enddo
     fret = 0 
  else
     write(LDT_logunit,*) 'ERR: Could not find file (I): ',trim(instfile)
     fret = -1
  endif
#endif
  
end subroutine retrieve_inst_ecmwfvars


!BOP
! !ROUTINE: retrieve_accum_ecmwfvars
! \label{retrieve_accum_ecmwfvars}
!
! !INTERFACE:
subroutine retrieve_accum_ecmwfvars(n, findex, avgfile1, avgfile2, glbdata1, &
                                    hr, fret)
  use LDT_coreMod,      only : LDT_rc, LDT_domain
  use LDT_logMod,       only : LDT_logunit, LDT_verify
  use ecmwf_forcingMod, only : ecmwf_struc
#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS:
  integer,   intent(in)        :: n 
  integer,   intent(in)        :: findex
  character(len=*), intent(in) :: avgfile1
  character(len=*), intent(in) :: avgfile2
  real                         :: glbdata1(9,LDT_rc%ngrid(n))
  integer                      :: hr
  integer, intent(inout)       :: fret
!
! !DESCRIPTION: 
! This routine opens the corresponding ECMWF data file to extract
! the specified variable, which represents an accumulated value.
! Should be used for shortwave, longwave, lsp, and cp.
! 
!EOP
  integer,parameter  :: N_ACCUM_VARS=4
  integer            :: necmwf
  integer            :: iret,gbret
  integer            :: c,r,iv, ftn,igrib
  integer            :: kk,nvars
  character(len=20)  :: shortName 
  real               :: missingValue 
  logical            :: file_exists
  character(len=4)   :: varname(N_ACCUM_VARS)
  real,  allocatable     :: f(:)
  logical*1, allocatable :: lb(:)
  integer            :: var_index
  logical            :: var_found
  logical            :: var_status(N_ACCUM_VARS)
  logical            :: pcp_flag
  real               :: temp_val
  real               :: glbdata2(9,LDT_rc%ngrid(n))
  integer            :: rel_index(N_ACCUM_VARS)
  real,dimension(LDT_rc%lnc(n), LDT_rc%lnr(n)) :: varfield
  logical            :: result1, result2  ! NaN check

  !=== set GRIB shortName 
  varname = (/ "ssrd","strd","lsp ","cp  " /)
  rel_index = (/ 3, 4, 8, 9 /) ! index of variable w.r.t. the
                               ! full list of 9 forcing variables

  necmwf = ecmwf_struc(n)%ncold*ecmwf_struc(n)%nrold

  varfield = 0
  var_status = .false.

#if (defined USE_GRIBAPI)
  !=== Set up to open file and retrieve specified field 
  fret = -1
  gbret = 0
  inquire(file=avgfile2,exist=file_exists) 
  if ( file_exists ) then 
     call grib_open_file(ftn,avgfile2,'r',iret)
     call LDT_verify(iret, 'failed to open file '//trim(avgfile2))

     call grib_count_in_file(ftn,nvars,iret)
     call LDT_verify(iret, 'error in grib_count_in_file in retrieve_avg_ecmwfvars')
     
     allocate(lb(necmwf))
     allocate(f(necmwf))

     do kk=1,nvars
        call grib_new_from_file(ftn,igrib,iret)
        call grib_get(igrib,'shortName',shortName)
        
        var_found = .false. 
        do iv=1,N_ACCUM_VARS
           if ( trim(shortName) == trim(varname(iv))) then
              var_found = .true. 
              var_index = rel_index(iv)
              var_status(iv) = .true. 
              exit
           endif
        enddo
        
        f = -9999.0
        call grib_get(igrib,'values',f,iret)
        call LDT_verify(iret, 'failed to get values in retrieve_accum_ecmwfvars')
        if ( iret .ne. 0 ) then  ! read not successful
         write(LDT_logunit,*) 'problem reading accum2 for ',shortName,iret
         call grib_release(igrib)
         deallocate(lb)
         deallocate(f)     
         fret = -1
         return 
        endif

        call grib_get(igrib,'missingValue',missingValue,iret)
        call LDT_verify(iret,'failed to get missingValue in retrieve_accum_ecmwfvars')

        call grib_release(igrib)
        call LDT_verify(iret,'error in grib_release in read_ecmwf')
        
        if(var_found) then
           lb = .false. 
           do c = 1, necmwf
              if ( f(c) .ne. missingValue ) then
                 lb(c) = .true.
              endif
           end do
           
           pcp_flag = .false.
           if(var_index.eq.8.or.var_index.eq.9) pcp_flag = .true.

           call interp_ecmwf(n,findex,pcp_flag,necmwf,f, &
                             lb,LDT_rc%gridDesc(n,:),    &
                             LDT_rc%lnc(n),LDT_rc%lnr(n),varfield)
           
           do r=1,LDT_rc%lnr(n)
              do c=1,LDT_rc%lnc(n)
                 if ( LDT_domain(n)%gindex(c,r) /= -1 ) then 
                    glbdata2(var_index,LDT_domain(n)%gindex(c,r)) = &
                       varfield(c,r)
                 endif
              enddo
           enddo
        endif

     enddo
     call grib_close_file(ftn)
  
     deallocate(lb)
     deallocate(f)     
     
     do kk=1,N_ACCUM_VARS
        if ( .not. var_status(kk) ) then 
           write(LDT_logunit,*) &
                'ERR: Could not retrieve all entries in file: ',trim(avgfile2)
           fret = -1
           return
        endif
     enddo
     fret = 0 
  else
     write(LDT_logunit,*) 'ERR: Could not find file (A2) : ',avgfile2   
     fret = -1
  endif
  gbret = gbret + fret

  var_status = .false.

  if ( .not. ((hr == 03) .or. (hr == 15)) ) then

     inquire(file=avgfile1,exist=file_exists) 
     if ( file_exists ) then 
        call grib_open_file(ftn,avgfile1,'r',iret)
        call LDT_verify(iret, 'failed to open file '//trim(avgfile1))
        
        call grib_count_in_file(ftn,nvars,iret)
        call LDT_verify(iret, 'error in grib_count_in_file in retrieve_accum_ecmwfvars')
     
        allocate(lb(necmwf))
        allocate(f(necmwf))

        do kk=1,nvars
           call grib_new_from_file(ftn,igrib,iret)
           call grib_get(igrib,'shortName',shortName)
           
           var_found = .false. 
           do iv=1,N_ACCUM_VARS
              if ( trim(shortName) == trim(varname(iv))) then
                 var_found = .true. 
                 var_index = rel_index(iv)
                 var_status(iv) = .true. 
                 exit
              endif
           enddo
           
           f = -9999.0
           call grib_get(igrib,'values',f,iret)
           call LDT_verify(iret, 'failed to get values in retrieve_accum_ecmwfvars')
           if ( iret .ne. 0 ) then   ! read not successful
            write(LDT_logunit,*) 'problem reading accum1 for ',shortName,iret
            call grib_release(igrib)
            deallocate(lb)
            deallocate(f)     
            fret = -1
            return 
           endif

           call grib_get(igrib,'missingValue',missingValue,iret)
           call LDT_verify(iret,'failed to get missingValue in retrieve_accum_ecmwfvars')
           call grib_release(igrib)
           call LDT_verify(iret,'error in grib_release in read_ecmwf')

           if(var_found) then
              lb = .false. 
              do c = 1, necmwf
                 if ( f(c) .ne. missingValue ) then
                    lb(c) = .true.
                 endif
              end do

              pcp_flag = .false.
              if(var_index.eq.8.or.var_index.eq.9) pcp_flag = .true.

              call interp_ecmwf(n,findex,pcp_flag,necmwf,f, &
                                lb,LDT_rc%gridDesc(n,:),    &
                                LDT_rc%lnc(n),LDT_rc%lnr(n),varfield)
              
              do r=1,LDT_rc%lnr(n)
                 do c=1,LDT_rc%lnc(n)
                    if ( LDT_domain(n)%gindex(c,r) /= -1 ) then 
                       glbdata1(var_index,LDT_domain(n)%gindex(c,r)) = &
                          varfield(c,r)
                    endif
                 enddo
              enddo
           endif   

        enddo
        call grib_close_file(ftn)
  
        deallocate(lb)
        deallocate(f)     
        
        do kk=1,N_ACCUM_VARS
           if ( .not. var_status(kk) ) then 
              write(LDT_logunit,*)                                  &
                   'ERR: Could not retrieve all entries in file: ', &
                   trim(avgfile1)
              fret = -1
              return
           endif
        enddo
        fret = 0 
     else
        write(LDT_logunit,*) 'ERR: Could not find file (A1): ',trim(avgfile1)
        fret = -1
     endif
     
  else
   do kk = 1, N_ACCUM_VARS
      iv = rel_index(kk)
      glbdata1(iv,:) = 0.0
   end do
  endif  ! not 03 or 15z
  gbret = gbret + fret
  
  if ( gbret .eq. 0 ) then  ! all fields valid
  do kk = 1, N_ACCUM_VARS
     iv = rel_index(kk)
     if ( iv == 3 ) then ! swd - convert accumulated field to rate
        if ((hr==3).or.(hr==15)) then
           glbdata1(iv,:) = glbdata2(iv,:) / (3.0*60*60)
        else if ((hr==6).or.(hr==18)) then
           glbdata1(iv,:) = glbdata2(iv,:)/ (6.0*60*60)
        else if ((hr==9).or.(hr==21)) then
           glbdata1(iv,:) = glbdata2(iv,:)/ (9.0*60*60)
        else if ((hr==0).or.(hr==12)) then
           glbdata1(iv,:) = glbdata2(iv,:) / (12.0*60*60)
        endif
     elseif ( iv == 8 .or. iv == 9 ) then ! lsp or cp
!=== Added because sometimes the f2 data is less than the f1 data
!=== creating negative precip. This ocurs only for small values of precip
!=== and the differences are on the order of 10E-4 so for intents and purposes 
!=== is not a substantial problem. This is an issue at the raw processing data 
!=== level issue that will need to be looked at more closely. Reprocessing 
!=== of entire dataset may be necessary if it is found that the problem can
!== be fixed. Will revisit. JG 3/10/2004.
!*** Added check and reset for NaN. H.B. 2/2/2015.
!*** Added check and reset for huge value. B.L. 2/5/2015.
        do c=1,LDT_rc%ngrid(n)
!           result1 = ISNAN(glbdata1(iv,c))
!           result2 = ISNAN(glbdata2(iv,c))
!           if ( result1 .eq. .true. ) glbdata1(iv,c) = 0.0
!           if ( result2 .eq. .true. ) glbdata2(iv,c) = 0.0
           temp_val = glbdata2(iv,c)-glbdata1(iv,c)
           if (temp_val>1000000.) temp_val = 0.0         
           if ( (glbdata2(iv,c)-glbdata1(iv,c)) < 0.0 ) then
              glbdata1(iv,c) = 0.0
           else
              glbdata1(iv,c) = (temp_val * 1000.0)/(3.0*60*60)
           endif
        enddo        
     else        
        glbdata1(iv,:) = (glbdata2(iv,:) - glbdata1(iv,:)) / (3.0*60*60)     
     endif
  enddo

  do c=1,LDT_rc%ngrid(n)
     glbdata1(8,c) = glbdata1(8,c) + glbdata1(9,c)
  enddo
  else
   fret = gbret
  endif  ! gbret == 0
#endif
  
end subroutine retrieve_accum_ecmwfvars

!BOP
! !ROUTINE: interp_ecmwf
! \label{interp_ecmwf}
!
! !INTERFACE:
subroutine interp_ecmwf(n, findex, pcp_flag, necmwf,f,lb,ldt_gds,nc,nr, &
     varfield)
! !USES:
  use LDT_coreMod,        only : LDT_rc, LDT_domain
  use ecmwf_forcingMod,   only : ecmwf_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
  logical, intent(in) :: pcp_flag
  integer, intent(in) :: necmwf
  real, intent(out)   :: f(necmwf)
  logical*1           :: lb(necmwf)
  real                :: ldt_gds(20)
  integer, intent(in) :: nc
  integer, intent(in) :: nr
  real, intent(out)   :: varfield(nc,nr)
!
! !DESCRIPTION:
!   This subroutine interpolates a given ECMWF field 
!   to the LDT grid. 
!  The arguments are: 
!  \begin{description}
! \item[n]
!  index of the nest
! \item[necmwf]
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
!-----------------------------------------------------------------------
! Setting interpolation options (ip=0,bilinear)
! (km=1, one parameter, ibi=1,use undefined bitmap
! (needed for soil moisture and temperature only)
! Use budget bilinear (ip=3) for precip forcing fields
!-----------------------------------------------------------------------
  mo = nc*nr
!-----------------------------------------------------------------------
! Initialize output bitmap. Important for soil moisture and temp.
!-----------------------------------------------------------------------
  lo = .true.
!-----------------------------------------------------------------------  
! Interpolate to LDT grid
!-----------------------------------------------------------------------  
!  if(LDT_rc%met_interp(findex).eq."bilinear") then
  if(LDT_rc%met_gridtransform(findex).eq."bilinear") then
     call bilinear_interp(ldt_gds,lb,f,lo,ldt1d,ecmwf_struc(n)%mi,mo,&
          LDT_domain(n)%lat, LDT_domain(n)%lon, &
          ecmwf_struc(n)%w111,ecmwf_struc(n)%w121,&
          ecmwf_struc(n)%w211,ecmwf_struc(n)%w221,&
          ecmwf_struc(n)%n111,ecmwf_struc(n)%n121,&
          ecmwf_struc(n)%n211,ecmwf_struc(n)%n221,LDT_rc%udef, iret)
!  elseif(LDT_rc%met_interp(findex).eq."budget-bilinear") then
  elseif(LDT_rc%met_gridtransform(findex).eq."budget-bilinear") then
     if (pcp_flag) then
        call conserv_interp(ldt_gds,lb,f,lo,ldt1d,ecmwf_struc(n)%mi,mo,&
             LDT_domain(n)%lat, LDT_domain(n)%lon, &
             ecmwf_struc(n)%w112,ecmwf_struc(n)%w122,&
             ecmwf_struc(n)%w212,ecmwf_struc(n)%w222,&
             ecmwf_struc(n)%n112,ecmwf_struc(n)%n122,&
             ecmwf_struc(n)%n212,ecmwf_struc(n)%n222,LDT_rc%udef, iret)
     else
        call bilinear_interp(ldt_gds,lb,f,lo,ldt1d,ecmwf_struc(n)%mi,mo,&
             LDT_domain(n)%lat, LDT_domain(n)%lon, &
             ecmwf_struc(n)%w111,ecmwf_struc(n)%w121,&
             ecmwf_struc(n)%w211,ecmwf_struc(n)%w221,&
             ecmwf_struc(n)%n111,ecmwf_struc(n)%n121,&
             ecmwf_struc(n)%n211,ecmwf_struc(n)%n221,LDT_rc%udef, iret)
     endif
  endif
!-----------------------------------------------------------------------    
! Create 2D array for main program. Also define a "soil" mask
! due to different geography between ECMWF & LDAS. For LDAS land 
! points not included in ECMWF geography dataset only.
!-----------------------------------------------------------------------    
  count1 = 0
  do j = 1, nr
     do i = 1, nc
        varfield(i,j) = ldt1d(i+count1)
     enddo
     count1 = count1 + nc
  enddo

end subroutine interp_ecmwf


