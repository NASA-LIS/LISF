!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: read_MODIS_RT_lai
! \label{read_MODIS_RT_lai}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!  20 Feb 2006: Sujay Kumar; Modified to support nesting
!
! !INTERFACE:
subroutine read_MODIS_RT_lai(n, forward_search, array, file_time_stamp)
! !USES:
  use ESMF
  use LIS_coreMod,        only : LIS_rc, LIS_domain, LIS_localPet
  use LIS_logMod,         only : LIS_logunit, LIS_getNextUnitNumber, &
       LIS_releaseUnitNumber, LIS_endrun
  use LIS_vegDataMod,    only : LIS_lai
  use LIS_fileIOMod,     only : LIS_readData
  use LIS_timeMgrMod, only: LIS_calendar
  use LIS_constantsMod,  only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 

  integer, intent(in)               :: n
  logical, intent(in)               :: forward_search
  type(ESMF_Time), intent (out)     :: file_time_stamp
  real, intent(inout)               :: array(LIS_rc%ntiles(n))

! !DESCRIPTION:
!  This subroutine retrieves the LAI data for the 
!  specified date and returns the values in the latlon projection
!  
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[mo]
!   time index (month or quarter)
!  \item[array]
!   output field with the retrieved LAI 
!  \end{description}
!
!EOP      
  integer*1, allocatable      :: tmparr(:,:)
  character(len=LIS_CONST_PATH_LEN) :: filename2
  integer :: yr, mo, dy 
  integer :: t, it, rc
  integer :: ftn, status, step
  character * 4 :: cyr
  character * 2 :: cmo, cdy 
  real        :: lai_gridDesc(6), lai
  logical     :: file_exists
  type(ESMF_Time)     :: time_lis, time
  type(ESMF_TimeInterval) :: deltaT

  allocate(tmparr(LIS_rc%lnc(n),LIS_rc%lnr(n)))
! MODIS RT LAI file example: 
!  input/LS_PARAMETERS/MODIS/MCD15A2.005/2004/lai.03.13.1gd1c
  call ESMF_TimeSet(time_lis, yy=LIS_rc%yr, &
       mm=LIS_rc%mo, &
       dd=LIS_rc%da, &
       h =LIS_rc%hr, &
       m =LIS_rc%mn, &
       s =LIS_rc%ss, &
       calendar = LIS_calendar, &
       rc = status)

  call ESMF_TimeGet(time_lis,yy=yr, mm=mo, dd=dy, calendar = LIS_calendar)

  !Handle forward or backward search (care with file for current day)
  if (forward_search.eqv..true.) then
     it=1  ! look beyond current time
     step=1
  else ! backward search
     it=0  ! start at current time & work backward
     step=-1
  endif
  
  do while (abs(it).le.8) ! search up to 8 days from current
     !  do it=1, 8  ! search forward for 8 days 
     call ESMF_TimeIntervalSet(deltaT,d=it,rc=rc)
     time = time_lis+deltaT
     call ESMF_TimeGet(time, yy=yr, mm=mo, dd=dy)
     
     write(cyr, '(I4.4)') yr
     write(cmo, '(I2.2)') mo
     write(cdy, '(I2.2)') dy
     filename2 = trim(LIS_lai(n)%laifile)//'/'//cyr//'/'//'lai.'//cmo//'.'//cdy//'.1gd1c' 
     
     inquire(file=trim(filename2), exist=file_exists)
     if(.not.file_exists) then 
        write(LIS_logunit,*) 'LAI map ',trim(filename2),' not found'
     else 
        exit
     end if
     it=it+step
  end Do
  
  if(.not.file_exists) then   ! did not find after 8 steps 
     write(LIS_logunit,*) 'LAI map ',trim(filename2),' not found'
     write(LIS_logunit,*) 'Program stopping ...'
     call LIS_endrun
  endif

  ftn = LIS_getNextUnitNumber()
  open(ftn, file=filename2, access='direct',status='old', &
       form="unformatted", recl=1)

  lai_gridDesc(1) = -59.875  !LIS_rc%gridDesc(n,4)
  lai_gridDesc(2) = -179.875 !LIS_rc%gridDesc(n,5)
  lai_gridDesc(3) = 89.875   !LIS_rc%gridDesc(n,7)
  lai_gridDesc(4) = 179.875  !LIS_rc%gridDesc(n,8)
  lai_gridDesc(5) = 0.25     !LIS_rc%gridDesc(n,9)
  lai_gridDesc(6) = 0.25     !LIS_rc%gridDesc(n,10)

  call read_lai_1gd1c(n,ftn,lai_gridDesc,tmparr)

  call LIS_releaseUnitNumber(ftn)

  do t=1,LIS_rc%ntiles(n)
     lai = tmparr(LIS_domain(n)%tile(t)%col,LIS_domain(n)%tile(t)%row) * 0.1
     ! reset undef to 0 
     if (lai .LT.0 .or. lai .GT. 10) lai = 0 
     array(t) = lai 
  enddo

  file_time_stamp = time
  deallocate(tmparr)

  write(LIS_logunit,*)'Read LAI File ',trim(filename2)
end subroutine read_MODIS_RT_lai

subroutine read_lai_1gd1c(n, ftn, gridDesc, char_array)
! !USES:
   use LIS_coreMod,  only : LIS_rc, LIS_domain, LIS_localPet
   !use gaussian_mod, only : gaussian_find_row, gaussian_find_col
   use map_utils,    only : ij_to_latlon
   use LIS_logMod,   only : LIS_abort, LIS_endrun, LIS_logunit

   implicit none
! !ARGUMENTS:
   integer,  INTENT(IN)           :: n
   integer,  INTENT(IN)           :: ftn
   integer*1, INTENT(INOUT)        :: char_array(LIS_rc%lnc(n),LIS_rc%lnr(n))
   real,     INTENT(IN)           :: gridDesc(6)
!
! !DESCRIPTION:
!   This routine retrieves a 2-d data from a binary direct access file.
!EOP
   real, pointer         :: rlat(:,:)
   real, pointer         :: rlon(:,:)
   integer*8    :: c, r
   integer*8    :: nc_dom
   integer*8    :: glnc, glnr
   integer*8    :: line1, line2, line
   integer      :: istat
   character*100  :: message  (20)
   character*1 :: temp_char

   allocate(rlat(LIS_rc%lnc(n),LIS_rc%lnr(n)))
   allocate(rlon(LIS_rc%lnc(n),LIS_rc%lnr(n)))

   if(trim(LIS_rc%lis_map_proj).eq."gaussian") then ! gaussian
#if 0
      line1 = gaussian_find_row(LIS_rc%gridDesc(n,4))  -   &
           gaussian_find_row(LIS_rc%gridDesc(n,44)) + 1
      line2 = gaussian_find_col(LIS_rc%gridDesc(n,5))  -   &
           gaussian_find_col(LIS_rc%gridDesc(n,45)) + 1

      do r=1,LIS_rc%lnr(n)
         do c=1,LIS_rc%lnc(n)
            glnc = line2+c-1
            glnr = line1+r-1
            line = (glnr-1)*nint(LIS_rc%gridDesc(n,42))+glnc
            read(ftn,rec=line,iostat=istat) char_array(c,r)
            if( istat .ne. 0 ) then
               message(1) = 'program:  LIS'
               message(2) = '  routine:  read2DData'
               message(3) = '  iostat != 0'
               call LIS_abort( message )
               call LIS_endrun
            endif
         enddo
      enddo
#else
      write(LIS_logunit, fmt=*) '[ERR] Gaussian support for MODIS RT LAI '// &
                                'has been disabled.'
#endif
   elseif(trim(LIS_rc%lis_map_proj).eq."latlon") then

      do r=1,LIS_rc%lnr(n)
         do c=1,LIS_rc%lnc(n)
            call ij_to_latlon(LIS_domain(n)%lisproj,float(c),float(r),&
                 rlat(c,r),rlon(c,r))
         enddo
      enddo

      nc_dom = nint((gridDesc(4)-gridDesc(2))/(gridDesc(5)))+1
      do r=1,LIS_rc%lnr(n)
         do c=1,LIS_rc%lnc(n)
            line1 = nint((rlat(c,r)-gridDesc(1))/gridDesc(6))+1
            line2 = nint((rlon(c,r)-gridDesc(2))/gridDesc(5))+1
            line = (line1-1)*nc_dom + line2
            read(ftn,rec=line) temp_char
            char_array(c,r)=ichar(temp_char)
         enddo
      enddo


   elseif(trim(LIS_rc%lis_map_proj).eq."UTM") then !utm
!rlat/rlon used here to store northing and easting
      do r=1,LIS_rc%lnr(n)
         do c=1,LIS_rc%lnc(n)
            rlat(c,r) = LIS_rc%gridDesc(n,4)+(r-1)*LIS_rc%gridDesc(n,9)
            rlon(c,r) = LIS_rc%gridDesc(n,5)+(c-1)*LIS_rc%gridDesc(n,9)
         enddo
      enddo

      nc_dom = gridDesc(4)

      do r=1,LIS_rc%lnr(n)
         do c=1,LIS_rc%lnc(n)
            line1 = nint((rlat(c,r)-gridDesc(2))/gridDesc(6))+1
            line2 = nint((rlon(c,r)-gridDesc(3))/gridDesc(6))+1
            line = (line1-1)*nc_dom +line2
           read(ftn,rec=line) char_array(c,r)
         enddo
      enddo
   else
      print*, 'This parameter projection is not supported...'
      print*, 'Program stopping ....'
      stop
   endif

   deallocate(rlat)
   deallocate(rlon)

end subroutine  read_lai_1gd1c



