!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: read_agrrad
! \label{read_agrrad}
!
! !REVISION HISTORY:
! 29Jul2005; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine read_agrrad(n,findex,order,yr,mo,da,hr)
! !USES: 
  use LIS_coreMod,        only : LIS_rc,LIS_domain
  use LIS_logMod,         only : LIS_logunit, LIS_endrun
  use LIS_timeMgrMod,     only : LIS_tick
  use LIS_metforcingMod,  only : LIS_forc
  use agrrad_forcingMod,  only : agrrad_struc

  implicit none
! !ARGUMENTS: 
  integer,   intent(in)    :: n 
  integer,   intent(in)    :: findex
  integer,   intent(in)    :: order
  integer,   intent(in)    :: yr,mo,da,hr
! 
! !DESCRIPTION: 
!  This routine opens the AGRMET files and reads the radiation 
!  fields. The fields are then spatially interpolated and 
!  gaps in the interpolation dues ot mismatches in landmask 
!  (between LIS and AGRMET) are filled. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]
!    index of the nest
!   \item[findex]
!    index of the forcing source
!   \item[order]
!    flag indicating which data to be read (order=1, read the previous 
!    hourly instance, order=2, read the next hourly instance)
!   \item[yr,mo,da,hr]
!    current year, month, day and hour. 
!  \end{description}
! 
!   The routines invoked are: 
!  \begin{description}
!   \item[LIS\_tick](\ref{LIS_tick}) \newline
!    determines the AGRMET data times
!   \item[agrrad\_filename](\ref{agrrad_filename}) \newline
!    generates the name of the AGRMET file to be read
!   \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!   \item[agrrad\_fillgaps](\ref{agrrad_fillgaps}) \newline
!    fills in the mismatches due to AGRMET and LIS masks
!   \item[retrieve\_agrrad\_variables](\ref{retrieve_agrrad_variables}) \newline
!    retrieves the agrrad variables 
!  \end{description}
!EOP

  integer,   parameter     :: nc = 720, nr=361
  character*100            :: agrradfile
  real                     :: var(nc,nr)
  real                     :: swd(nc*nr)
  logical*1                :: lb(nc*nr)
  real                     :: lwd(nc*nr)
  integer                  :: iret
  logical                  :: exists1
  integer                  :: c,r
  integer                  :: swd_index,lwd_index
  logical*1                :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  real                     :: go(LIS_rc%lnc(n)*LIS_rc%lnr(n)) 
  real                     :: varfield(LIS_rc%lnc(n),LIS_rc%lnr(n))
  integer                  :: gindex
  integer                  :: yr1, mo1,da1,hr1,mn1,ss1,doy1
  real                     :: ts1, gmt1
  real*8                   :: timenow
  integer                  :: error_swd, error_lwd

  yr1 = LIS_rc%yr  !current time
  mo1 = LIS_rc%mo
  da1 = LIS_rc%da
  hr1 = LIS_rc%hr
  mn1 = LIS_rc%mn
  ss1 = 0
  ts1 = 0
  call LIS_tick( timenow, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )
  write(LIS_logunit,*) "MSG: read_agrrad: timenow = ", timenow
  write(LIS_logunit,*) "MSG: read_agrrad: changetime1 = ", agrrad_struc(n)%changetime1
  write(LIS_logunit,*) "MSG: read_agrrad: changetime2 = ", agrrad_struc(n)%changetime2
  if(timenow.lt.agrrad_struc(n)%changetime1) then 
     swd_index = 12
     lwd_index = 13
  elseif((timenow.ge.agrrad_struc(n)%changetime1).and.&
       (timenow.lt.agrrad_struc(n)%changetime2)) then 
     swd_index = 30
     lwd_index = 31
  else
     swd_index = 31
     lwd_index = 32
  endif
  write(LIS_logunit,*) "MSG: read_agrrad: swd_index = ", swd_index
  write(LIS_logunit,*) "MSG: read_agrrad: lwd_index = ", lwd_index
  call agrrad_filename(agrradfile,agrrad_struc(n)%agrdir,&
       yr,mo,da,hr)
  
  error_swd = -1 ; error_lwd = -1 ! assume failure

  lb = .true.
  inquire(file=trim(agrradfile),exist=exists1)

  if(exists1) then 
     write(LIS_logunit,*) 'Reading AGRMET file ',trim(agrradfile)

     call retrieve_agrrad_variables(agrradfile, error_swd, swd_index, var)
     do r=1,nr
        do c=1,nc
           if(c.ge.1.and.c.le.360) then
              swd(c+(r-1)*nc) = var(c+360,361-r+1)
           else
              swd(c+(r-1)*nc) = var(c-360,361-r+1)
           endif
           if(swd(c+(r-1)*nc).eq.-1E30) lb(c+(r-1)*nc) = .false.
        enddo
     enddo

     call retrieve_agrrad_variables(agrradfile, error_lwd, lwd_index, var)
     do r=1,nr
        do c=1,nc
           if(c.ge.1.and.c.le.360) then
              lwd(c+(r-1)*nc) = var(c+360,361-r+1)
           else
              lwd(c+(r-1)*nc) = var(c-360,361-r+1)
           endif
        enddo
     enddo
  endif

  if ( error_swd == 0 .and. error_lwd == 0 )  then 
       
     if(LIS_rc%met_interp(findex).eq."bilinear") then 
        call bilinear_interp(LIS_rc%gridDesc(n,:),lb,swd,lo,go,&
             nc*nr,LIS_rc%lnc(n)*LIS_rc%lnr(n),&
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             agrrad_struc(n)%w111,agrrad_struc(n)%w121,&
             agrrad_struc(n)%w211,agrrad_struc(n)%w221,&
             agrrad_struc(n)%n111,agrrad_struc(n)%n121,&
             agrrad_struc(n)%n211,agrrad_struc(n)%n221,&
             LIS_rc%udef,iret)
     endif
     
     varfield = -9999.0
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           varfield(c,r) = go(c+(r-1)*LIS_rc%lnc(n))
        enddo
     enddo
     call agrrad_fillgaps(n,varfield)
     
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           gindex = LIS_domain(n)%gindex(c,r)
           if(gindex.ne.-1) then 
              if ( order == 1 ) then
                 agrrad_struc(n)%metdata1(1,gindex) = varfield(c,r)
              else
                 agrrad_struc(n)%metdata2(1,gindex) = varfield(c,r)
              endif
           endif
        enddo
     enddo
     

     if(LIS_rc%met_interp(findex).eq."bilinear") then 
        call bilinear_interp(LIS_rc%gridDesc(n,:),lb,lwd,lo,go,&
             nc*nr,LIS_rc%lnc(n)*LIS_rc%lnr(n),&
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             agrrad_struc(n)%w111,agrrad_struc(n)%w121,&
             agrrad_struc(n)%w211,agrrad_struc(n)%w221,&
             agrrad_struc(n)%n111,agrrad_struc(n)%n121,&
             agrrad_struc(n)%n211,agrrad_struc(n)%n221,&
             LIS_rc%udef,iret)
     endif
     
     varfield = -9999.0
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           varfield(c,r) = go(c+(r-1)*LIS_rc%lnc(n))
        enddo
     enddo
     
     call agrrad_fillgaps(n,varfield)
     
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           gindex = LIS_domain(n)%gindex(c,r)
           if(gindex.ne.-1) then 
              if ( order == 1 ) then
                 agrrad_struc(n)%metdata1(2,gindex) = varfield(c,r)
              else
                 agrrad_struc(n)%metdata2(2,gindex) = varfield(c,r)
              endif
           endif
        enddo
     enddo
  else
     if ( order == 1 ) then
        agrrad_struc(n)%metdata1 = LIS_rc%udef
     else
        agrrad_struc(n)%metdata2 = LIS_rc%udef
     endif
  endif

 end subroutine read_agrrad


!BOP
! 
! !ROUTINE: agrrad_filename
! \label{agrrad_filename}
! 
! !INTERFACE: 
 subroutine agrrad_filename(name,rootdir,yr,mo,da,hr)

   implicit none
! !ARGUMENTS:   
  
   character(len=*)          :: name
   character(len=*)          :: rootdir
   integer, intent(IN)       :: yr
   integer, intent(IN)       :: mo
   integer, intent(IN)       :: da
   integer, intent(IN)       :: hr
! 
! !DESCRIPTION: 
!  This routines generates the name of the surface temperature file
!  by appending the hemisphere and timestamps to the root directory. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[name] filename to be generated
!   \item[rootdir] name of the root directory that contains the forcing
!   \item[yr] year of data
!   \item[mo] month of  data
!   \item[da] day of data
!   \item[hr] hour of data
!  \end{description}
!
!EOP
  
   character(len=6)    :: ftime1
   character(len=10)   :: ftime2
   
   write(unit=ftime1,fmt='(i4,i2.2)') yr,mo
   write(unit=ftime2,fmt='(i4,i2.2,i2.2,i2.2)') yr,mo,da,hr
   
   name = trim(rootdir)//'/'//ftime1//'/agrmet.grib.03hr.'//ftime2
 end subroutine agrrad_filename
 

!BOP
! !ROUTINE: agrrad_fillgaps
!  \label{agrrad_fillgaps}
! 
! !INTERFACE:
subroutine agrrad_fillgaps(n,varfield)
! !USES:
  use LIS_coreMod,       only : LIS_rc,LIS_domain
  use LIS_logMod,      only   : LIS_logunit, LIS_endrun
  use agrrad_forcingMod, only : agrrad_struc

  implicit none
! !USES: 
  integer, intent(in)    :: n
  real, intent(inout)    :: varfield(LIS_rc%lnc(n),LIS_rc%lnr(n)) 
! !DESCRIPTION:
!   This subroutine fills in invalid grid points introduced due to 
!   reprojection from PS to lat/lon. This routine assumes that the undef
!   or invalid value is the LIS undefined value. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[ip]
!    interpolation option
!  \item[varfield]
!    updated output field
!  \end{description}
!
!EOP
  integer                :: c,r
  logical                :: foundPt
  integer                :: i,j,str,enr,stc,enc,kk
  integer                :: try

  try = 0 
  if(agrrad_struc(n)%fillflag1) then !This will be done once 
     agrrad_struc(n)%smask1 = 0
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if((LIS_domain(n)%gindex(c,r).ne.-1).and.&
                varfield(c,r).eq.LIS_rc%udef) then !mismatch
              agrrad_struc(n)%smask1(c,r) = 1
           endif
        enddo
     enddo
     agrrad_struc(n)%fillflag1 = .false. 
  endif
     
  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(agrrad_struc(n)%smask1(c,r).eq.1) then 
           foundPt = .false.
           kk = 1
           try = 0 
           do while(.not.foundPt) 
              try = try +1
              str = max(r-kk,1)
              enr = min(r+kk,LIS_rc%lnr(n))
              stc = max(c-kk,1)
              enc = min(c+kk,LIS_rc%lnc(n))
              do j=str,enr
                 do i=stc,enc
                    if(LIS_domain(n)%gindex(i,j).ne.-1&
                         .and.agrrad_struc(n)%smask1(i,j).ne.1) then 
                       varfield(c,r) = varfield(i,j)
                       foundPt = .true.
                       exit
                    endif
                 enddo
              enddo
              kk = kk+1
              if(try.gt.100) then 
                 write(LIS_logunit,*) 'AGRMET fillgaps failed, stopping..',try,kk,c,r
                 call LIS_endrun()
              endif
           enddo
        endif
     enddo
  enddo
     
end subroutine agrrad_fillgaps

!BOP
! !ROUTINE: retrieve_agrrad_variables
! \label{retrieve_agrrad_variables}
!
! !REVISION HISTORY:
!  04 Dec 2013: James Geiger: initial implementation
!
! !INTERFACE:
subroutine retrieve_agrrad_variables(fname, error, vindex, var)
! !USES:  
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only : LIS_logunit, LIS_endrun, LIS_verify

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
  integer, parameter             :: nc=720, nr=361
! !ARGUMENTS:
  character(len=100), intent(in) :: fname
  integer, intent(out)           :: error
  integer, intent(in)            :: vindex
  real, intent(out)              :: var(nc,nr)
!
! !DESCRIPTION:
!  This routine retrieves the AGRRAD forcing variable specified by
!  vindex from the file specified by fname.
!
!  The arguments are: 
!  \begin{description}
!  \item[fname]
!   name of the AGRRAD forcing file
!  \item[error]
!   return status of this subroutine
!   error = 0 indicates success
!   error = -1 indicates failure
!  \item[vindex]
!   index of the variable to read
!  \item[var]
!   array to contain AGRRAD forcing variable read from file
!  \end{description}
! 
!EOP

!==== Local Variables=======================
  
#if (defined USE_GRIBAPI)
  integer :: ftn, igrib
  integer :: i, iret
  real    :: missingValue
  real    :: var1(nc*nr)

  error = 0 ! assume a successful read

  call grib_open_file(ftn,trim(fname),'r',iret)
  if ( iret /= 0 ) then
     write(LIS_logunit,*) 'ERR: retrieve_agrrad_variables: '// &
                          'Could not open '                 // &
                          trim(fname)

     error = -1
     return
  endif

  do i = 1, vindex-1
     call grib_new_from_file(ftn, igrib, iret)
     call grib_release(igrib,iret)
  enddo

  call grib_new_from_file(ftn, igrib, iret)

  var  = -9999.0
  var1 = -9999.0
  call grib_get(igrib,'values',var1,iret)

  if(iret /= 0) then 
     write(LIS_logunit,*) 'ERR: retrive_agrrad_variables: ' //     &
                          'Could not retrieve values in file: ' // &
                          trim(fname)
     write(LIS_logunit,*) 'ERR: retrive_agrrad_variables: ' //       &
                          'Could not retrieve values for forcing: ', &
                          vindex
     error = -1
     return
  endif

  var = reshape(var1, (/nc,nr/))

  call grib_get(igrib,'missingValue',missingValue,iret)
  if ( iret /= 0 ) then
     !write(LIS_logunit,*) 'ERR: retrieve_agrrad_variables: ' // &
     !                     'Could not get missingValue'
     missingValue = LIS_rc%udef
  endif

  call grib_release(igrib,iret)
  call LIS_verify(iret, 'ERR: retrieve_agrrad_variabes: ' // &
                        'Could not release igrib')

  call grib_close_file(ftn)
#else
  write(LIS_logunit,*) 'ERR: retrive_agrrad_variables: ' //              &
                       'AGRRAD support requires the GRIB_API library. ', &
                       'Please recompile LIS.'

  call LIS_endrun
#endif
end subroutine retrieve_agrrad_variables
