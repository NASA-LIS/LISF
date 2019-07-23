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
! !ROUTINE: read_pildas
! \label{read_pildas}
! 
! !REVISION HISTORY:
! 25 Apr 2013: Sujay Kumar, Initial code
! 14 Jul 2016: Mahdi Navari - Modified for PILDAS
! !INTERFACE:      
subroutine read_pildas(n, findex, order, name,ferror)
! !USES:
  use LIS_coreMod, only        : LIS_rc, LIS_domain
  use LIS_logMod, only         : LIS_logunit, LIS_verify, LIS_warning
  use LIS_metforcingMod, only  : LIS_forc
  use pildas_forcingMod,only : pildas_struc
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  implicit none
! !ARGUMENTS: 
  integer, intent(in)      :: findex
  integer, intent(in)      :: n
  integer, intent(in)      :: order
  character*80, intent(in) :: name
  integer, intent(out)     :: ferror
!
! !DESCRIPTION:
!  For the given time, reads the forcing data from the 
!  GEOS5 file, transforms into LIS forcing 
!  parameters and interpolates to the LIS domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    index of the forcing
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous 
!    hourly instance, order=2, read the next hourly instance)
!  \item[name]
!    name of the file to be read
!  \item[ferror]
!    return error flag (0-fail, 1-success)
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_pildas](\ref{interp_pildas}) \newline
!    Performs spatial interpolation of GEOS5 forecast data to the LIS grid
!  \end{description}

!EOP
  integer             :: ftn
  logical             :: file_exists  
  integer             :: c,r
  integer             :: t2mId,q2mId,swdnId,lwdnId,u10mId,v10mId
  integer             :: psurfId,prectotId,precconId,precsnoId
  real                :: t2m(pildas_struc(n)%nc,pildas_struc(n)%nr)
  real                :: q2m(pildas_struc(n)%nc,pildas_struc(n)%nr)
  real                :: swdn(pildas_struc(n)%nc,pildas_struc(n)%nr)
  real                :: lwdn(pildas_struc(n)%nc,pildas_struc(n)%nr)
  real                :: u10m(pildas_struc(n)%nc,pildas_struc(n)%nr)
  real                :: v10m(pildas_struc(n)%nc,pildas_struc(n)%nr)
  real                :: psurf(pildas_struc(n)%nc,pildas_struc(n)%nr)
  real                :: prectot(pildas_struc(n)%nc,pildas_struc(n)%nr)
  real                :: preccon(pildas_struc(n)%nc,pildas_struc(n)%nr)
  real                :: precsno(pildas_struc(n)%nc,pildas_struc(n)%nr)
  real                :: varfield(LIS_rc%lnc(n),LIS_rc%lnr(n))

  ferror = 1
  
  varfield = 0 
  ferror = 1

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  inquire (file=name, exist=file_exists)
  if (file_exists) then      

     call LIS_verify(nf90_open(path=trim(name),mode=NF90_NOWRITE,ncid=ftn),&
          'nf90_open failed in read_pildas')

     if(pildas_struc(n)%uselml.eq.1) then 
        call LIS_verify(nf90_inq_varid(ftn,'TLML',t2mId),&
             'nf90_inq_varid failed for TLML in read_pildas')
        call LIS_verify(nf90_inq_varid(ftn,'QLML',q2mId),&
             'nf90_inq_varid failed for QLML in read_pildas')
        call LIS_verify(nf90_inq_varid(ftn,'ULML',u10mId),&
             'nf90_inq_varid failed for ULML in read_pildas')
        call LIS_verify(nf90_inq_varid(ftn,'VLML',v10mId),&
             'nf90_inq_varid failed for VLML in read_pildas')
     else
        call LIS_verify(nf90_inq_varid(ftn,'T2M',t2mId),&
             'nf90_inq_varid failed for T2M in read_pildas')
        call LIS_verify(nf90_inq_varid(ftn,'Q2M',q2mId),&
             'nf90_inq_varid failed for Q2M in read_pildas')
        call LIS_verify(nf90_inq_varid(ftn,'U10M',u10mId),&
             'nf90_inq_varid failed for U10M in read_pildas')
        call LIS_verify(nf90_inq_varid(ftn,'V10M',v10mId),&
             'nf90_inq_varid failed for V10M in read_pildas')
     endif

     call LIS_verify(nf90_inq_varid(ftn,'SWDN',swdnId),&
             'nf90_inq_varid failed for SWDN in read_pildas')
     call LIS_verify(nf90_inq_varid(ftn,'LWDN',lwdnId),&
          'nf90_inq_varid failed for LWDN in read_pildas')
     call LIS_verify(nf90_inq_varid(ftn,'PSURF',psurfId),&
          'nf90_inq_varid failed for PSURF in read_pildas')
     call LIS_verify(nf90_inq_varid(ftn,'PRECTOT',prectotId),&
          'nf90_inq_varid failed for PRECTOT in read_pildas')
     call LIS_verify(nf90_inq_varid(ftn,'PRECCON',precconId),&
          'nf90_inq_varid failed for PRECCON in read_pildas')
     call LIS_verify(nf90_inq_varid(ftn,'PRECSNO',precsnoId),&
          'nf90_inq_varid failed for PRECSNO in read_pildas')
     
     call LIS_verify(nf90_get_var(ftn,t2mid,t2m),&
          'nf90_get_var failed for T2M in read_pildas')
     call LIS_verify(nf90_get_var(ftn,q2mid,q2m),&
          'nf90_get_var failed for Q2M in read_pildas')
     call LIS_verify(nf90_get_var(ftn,swdnid,swdn),&
          'nf90_get_var failed for SWDN in read_pildas')     
     call LIS_verify(nf90_get_var(ftn,lwdnid,lwdn),&
          'nf90_get_var failed for LWDN in read_pildas')
     call LIS_verify(nf90_get_var(ftn,u10mid,u10m),&
          'nf90_get_var failed for U10M in read_pildas')
     call LIS_verify(nf90_get_var(ftn,v10mid,v10m),&
          'nf90_get_var failed for V10M in read_pildas')
     call LIS_verify(nf90_get_var(ftn,psurfid,psurf),&
          'nf90_get_var failed for PSURF in read_pildas')
     call LIS_verify(nf90_get_var(ftn,prectotid,prectot),&
          'nf90_get_var failed for PRECTOT in read_pildas')
     call LIS_verify(nf90_get_var(ftn,precconid,preccon),&
          'nf90_get_var failed for PRECCON in read_pildas')
     call LIS_verify(nf90_get_var(ftn,precsnoid,precsno),&
          'nf90_get_var failed for PRECSNO in read_pildas')

     call LIS_verify(nf90_close(ftn))

     call interp_pildas(n, findex,t2m, varfield)
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 pildas_struc(n)%metdata1(1,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 pildas_struc(n)%metdata2(1,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo

     call interp_pildas(n, findex,q2m, varfield)
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 pildas_struc(n)%metdata1(2,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 pildas_struc(n)%metdata2(2,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo

     call interp_pildas(n, findex,swdn, varfield)
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 pildas_struc(n)%metdata1(3,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 pildas_struc(n)%metdata2(3,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo

     call interp_pildas(n, findex,lwdn, varfield)
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 pildas_struc(n)%metdata1(4,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 pildas_struc(n)%metdata2(4,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo

     call interp_pildas(n, findex,u10m, varfield)

     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 pildas_struc(n)%metdata1(5,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 pildas_struc(n)%metdata2(5,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo

     call interp_pildas(n, findex,v10m, varfield)

     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 pildas_struc(n)%metdata1(6,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 pildas_struc(n)%metdata2(6,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo


     call interp_pildas(n, findex,psurf, varfield)
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 pildas_struc(n)%metdata1(7,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 pildas_struc(n)%metdata2(7,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo

     call interp_pildas(n, findex,prectot, varfield)
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 pildas_struc(n)%metdata1(8,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 pildas_struc(n)%metdata2(8,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo

     call interp_pildas(n, findex,preccon, varfield)
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 pildas_struc(n)%metdata1(9,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 pildas_struc(n)%metdata2(9,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo

     call interp_pildas(n, findex,precsno, varfield)
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 pildas_struc(n)%metdata1(10,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 pildas_struc(n)%metdata2(10,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo

  else
     write(LIS_logunit,*) &
          'Could not find file: ',trim(name)
     ferror = 0
  endif
#endif     

end subroutine read_pildas


