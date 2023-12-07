!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_geis
!  \label{read_geis}
!
! !REVISION HISTORY:
! 15 Nov 2023; Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine read_geis(n, kk, findex, order, month, name,ferror)
! !USES:
  use LIS_coreMod
  use LIS_logMod, only         : LIS_logunit, LIS_verify, LIS_warning
  use LIS_metforcingMod, only  : LIS_forc
  use geis_forcingMod, only  : geis_struc
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in)       :: n
  integer, intent(in)       :: kk      ! Forecast member index
  integer, intent(in)       :: findex  ! Forcing index
  integer, intent(in)       :: order
  integer, intent(out)      :: month
  character(len=*), intent(in) :: name
  integer, intent(out)      :: ferror
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  the Global EIS data, transforms into 11 LIS forcing 
!  parameters 
!
!  The arguments are: 
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous 
!    hourly instance, order=2, read the next hourly instance)
!  \item[n]
!    index of the nest
!  \item[name]
!    name of the hourly GEIS forecast file
!  \item[ferror]
!    flag to indicate success of the call (=1 indicates success)
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_geis](\ref{interp_geis}) \newline
!    spatially interpolates a GEIS variable
!  \end{description}
!EOP

  integer                   :: ftn
  integer                   :: geis,varid
  integer                   :: k,t,c,r,iret,rc
  logical                   :: file_exists
  real                      :: varfield(LIS_rc%lnc(n),LIS_rc%lnr(n))

  ferror = 1

#if (defined USE_NETCDF4)
  geis = (geis_struc(n)%ncold*geis_struc(n)%nrold)
  
  varfield = 0 
  ferror = 1

  inquire (file=name, exist=file_exists)
  if (file_exists) then      


     call LIS_verify(nf90_open(path=trim(name), mode=NF90_NOWRITE, &
          ncid=ftn), 'nf90_open failed in read_geis')
          
     call LIS_verify(nf90_inq_varid(ftn,'Tair',varId), &
          'nf90_inq_varid failed for Tair in read_geis')
     call LIS_verify(nf90_get_var(ftn,varId, varfield,&
          start=(/geis_struc(n)%c_off,&
          geis_struc(n)%r_off/),&
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/)),&
          'nf90_get_var failed for Tair in read_geis')             

     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 geis_struc(n)%metdata1(kk,1,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 geis_struc(n)%metdata2(kk,1,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo

     call LIS_verify(nf90_inq_varid(ftn,'Qair',varId), &
          'nf90_inq_varid failed for Qair in read_geis')
     call LIS_verify(nf90_get_var(ftn,varId, varfield,&
          start=(/geis_struc(n)%c_off,&
          geis_struc(n)%r_off/),&
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/)),&
          'nf90_get_var failed for Qair in read_geis')             

     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 geis_struc(n)%metdata1(kk,2,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 geis_struc(n)%metdata2(kk,2,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo


     call LIS_verify(nf90_inq_varid(ftn,'SWdown',varId), &
          'nf90_inq_varid failed for SWdown in read_geis')
     call LIS_verify(nf90_get_var(ftn,varId, varfield,&
          start=(/geis_struc(n)%c_off,&
          geis_struc(n)%r_off/),&
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/)),&
          'nf90_get_var failed for SWdown in read_geis')             

     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 geis_struc(n)%metdata1(kk,3,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 geis_struc(n)%metdata2(kk,3,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo


     call LIS_verify(nf90_inq_varid(ftn,'LWdown',varId), &
          'nf90_inq_varid failed for LWdown in read_geis')
     call LIS_verify(nf90_get_var(ftn,varId, varfield,&
          start=(/geis_struc(n)%c_off,&
          geis_struc(n)%r_off/),&
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/)),&
          'nf90_get_var failed for LWdown in read_geis')             

     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 geis_struc(n)%metdata1(kk,4,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 geis_struc(n)%metdata2(kk,4,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo
          

     call LIS_verify(nf90_inq_varid(ftn,'Wind',varId), &
          'nf90_inq_varid failed for Wind in read_geis')
     call LIS_verify(nf90_get_var(ftn,varId, varfield,&
          start=(/geis_struc(n)%c_off,&
          geis_struc(n)%r_off/),&
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/)),&
          'nf90_get_var failed for Wind in read_geis')             

     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 geis_struc(n)%metdata1(kk,5,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 geis_struc(n)%metdata2(kk,5,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo

     call LIS_verify(nf90_inq_varid(ftn,'Psurf',varId), &
          'nf90_inq_varid failed for Psurf in read_geis')
     call LIS_verify(nf90_get_var(ftn,varId, varfield,&
          start=(/geis_struc(n)%c_off,&
          geis_struc(n)%r_off/),&
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/)),&
          'nf90_get_var failed for Psurf in read_geis')             

     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 geis_struc(n)%metdata1(kk,6,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 geis_struc(n)%metdata2(kk,6,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo

     call LIS_verify(nf90_inq_varid(ftn,'Rainf',varId), &
          'nf90_inq_varid failed for Rainf in read_geis')
     call LIS_verify(nf90_get_var(ftn,varId, varfield,&
          start=(/geis_struc(n)%c_off,&
          geis_struc(n)%r_off/),&
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/)),&
          'nf90_get_var failed for Rainf in read_geis')             

     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 geis_struc(n)%metdata1(kk,7,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 geis_struc(n)%metdata2(kk,7,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo
      
     
  else
     write(LIS_logunit,*) &
          '[ERR] Could not find file: ',trim(name)
     ferror = 0
  endif
#endif     

end subroutine read_geis

