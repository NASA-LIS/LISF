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
! !ROUTINE: read_geos5fcst
!  \label{read_geos5fcst}
!
! !REVISION HISTORY:
! 7 Mar 2013: Sujay Kumar, initial specification
!
! !INTERFACE:
subroutine read_geos5fcst(n, m, findex, order, name, ferror)
! !USES:
  use LIS_coreMod,          only : LIS_rc, LIS_domain
  use LIS_logMod,           only : LIS_logunit, LIS_verify, LIS_warning
  use LIS_metforcingMod,    only : LIS_forc
  use geos5fcst_forcingMod, only : geos5fcst_struc
  use LIS_constantsMod,     only : LIS_CONST_PATH_LEN
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  implicit none
! !ARGUMENTS: 
  integer, intent(in)       :: findex
  integer, intent(in)       :: n
  integer, intent(in)       :: m
  integer, intent(in)       :: order
  character(len=LIS_CONST_PATH_LEN), intent(in) :: name
  integer, intent(out)      :: ferror
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
!  \item[interp\_geos5fcst](\ref{interp_geos5fcst}) \newline
!    Performs spatial interpolation of GEOS5 forecast data to the LIS grid
!  \end{description}

!EOP
  integer             :: ftn
  logical             :: file_exists  
  integer             :: c,r
  integer             :: tlmlId,qlmlId,swgdnId,lwgabId,speedlmlId,psId,prectotId
  integer             :: precconId,precsnoId,swlandId
  !SB
  integer             :: hlmlId, pardrId, pardfId
  !EB
  real                :: tlml(geos5fcst_struc(n)%nc,geos5fcst_struc(n)%nr)
  real                :: qlml(geos5fcst_struc(n)%nc,geos5fcst_struc(n)%nr)
  real                :: swgdn(geos5fcst_struc(n)%nc,geos5fcst_struc(n)%nr)
  real                :: lwgab(geos5fcst_struc(n)%nc,geos5fcst_struc(n)%nr)
  real                :: speedlml(geos5fcst_struc(n)%nc,geos5fcst_struc(n)%nr)
  real                :: ps(geos5fcst_struc(n)%nc,geos5fcst_struc(n)%nr)
  real                :: prectot(geos5fcst_struc(n)%nc,geos5fcst_struc(n)%nr)
  real                :: precsno(geos5fcst_struc(n)%nc,geos5fcst_struc(n)%nr)
  real                :: preccon(geos5fcst_struc(n)%nc,geos5fcst_struc(n)%nr)
  !SB
  real                :: hlml(geos5fcst_struc(n)%nc,geos5fcst_struc(n)%nr)
  real                :: pardr(geos5fcst_struc(n)%nc,geos5fcst_struc(n)%nr)
  real                :: pardf(geos5fcst_struc(n)%nc,geos5fcst_struc(n)%nr)
  real                :: swland(geos5fcst_struc(n)%nc,geos5fcst_struc(n)%nr)
  !EB
  real                :: varfield(LIS_rc%lnc(n),LIS_rc%lnr(n))

  ferror = 1
  
  varfield = 0 
  ferror = 1

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  inquire (file=name, exist=file_exists)
  if (file_exists) then      

     call LIS_verify(nf90_open(path=trim(name),mode=NF90_NOWRITE,ncid=ftn),&
          'nf90_open failed in read_geos5fcst')

     call LIS_verify(nf90_inq_varid(ftn,'TLML',tlmlId),&
          'nf90_inq_varid failed for TLML in read_geos5fcst')
     call LIS_verify(nf90_inq_varid(ftn,'QLML',qlmlId),&
          'nf90_inq_varid failed for QLML in read_geos5fcst')
     call LIS_verify(nf90_inq_varid(ftn,'SWGDN',swgdnId),&
          'nf90_inq_varid failed for SWGDN in read_geos5fcst')
     call LIS_verify(nf90_inq_varid(ftn,'LWGAB',lwgabId),&
          'nf90_inq_varid failed for LWGAB in read_geos5fcst')
     call LIS_verify(nf90_inq_varid(ftn,'SPEEDLML',speedlmlId),&
          'nf90_inq_varid failed for SPEEDLML in read_geos5fcst')
     call LIS_verify(nf90_inq_varid(ftn,'PS',psId),&
          'nf90_inq_varid failed for PS in read_geos5fcst')
     call LIS_verify(nf90_inq_varid(ftn,'PRECTOT',prectotId),&
          'nf90_inq_varid failed for PRECTOT in read_geos5fcst')
     call LIS_verify(nf90_inq_varid(ftn,'PRECSNO',precsnoId),&
          'nf90_inq_varid failed for PRECSNO in read_geos5fcst')
     call LIS_verify(nf90_inq_varid(ftn,'PRECCON',precconId),&
          'nf90_inq_varid failed for PRECCON in read_geos5fcst')
!SB
     call LIS_verify(nf90_inq_varid(ftn,'HLML',hlmlId),&
          'nf90_inq_varid failed for HLML in read_geos5fcst')
     call LIS_verify(nf90_inq_varid(ftn,'PARDR',pardrId),&
          'nf90_inq_varid failed for PARDR in read_geos5fcst')
     call LIS_verify(nf90_inq_varid(ftn,'PARDF',pardfId),&
          'nf90_inq_varid failed for PARDF in read_geos5fcst')
     call LIS_verify(nf90_inq_varid(ftn,'SWLAND',swlandId),&
          'nf90_inq_varid failed for SWLAND in read_geos5fcst')
!EB

     call LIS_verify(nf90_get_var(ftn,tlmlid,tlml),&
          'nf90_get_var failed for TLML in read_geos5fcst')
     call LIS_verify(nf90_get_var(ftn,qlmlid,qlml),&
          'nf90_get_var failed for QLML in read_geos5fcst')
     call LIS_verify(nf90_get_var(ftn,swgdnid,swgdn),&
          'nf90_get_var failed for SWGDN in read_geos5fcst')     
     call LIS_verify(nf90_get_var(ftn,lwgabid,lwgab),&
          'nf90_get_var failed for LWGAB in read_geos5fcst')
     call LIS_verify(nf90_get_var(ftn,speedlmlid,speedlml),&
          'nf90_get_var failed for SPEEDLML in read_geos5fcst')
     call LIS_verify(nf90_get_var(ftn,psid,ps),&
          'nf90_get_var failed for PS in read_geos5fcst')
     call LIS_verify(nf90_get_var(ftn,prectotid,prectot),&
          'nf90_get_var failed for PRECTOT in read_geos5fcst')
     call LIS_verify(nf90_get_var(ftn,precsnoid,precsno),&
          'nf90_get_var failed for PRECSNO in read_geos5fcst')
     call LIS_verify(nf90_get_var(ftn,precconid,preccon),&
          'nf90_get_var failed for PRECCON in read_geos5fcst')
!SB
     call LIS_verify(nf90_get_var(ftn,hlmlid,hlml),&
          'nf90_get_var failed for HLML in read_geos5fcst')
     call LIS_verify(nf90_get_var(ftn,pardrid,pardr),&
          'nf90_get_var failed for PARDR in read_geos5fcst')
     call LIS_verify(nf90_get_var(ftn,pardfid,pardf),&
          'nf90_get_var failed for PARDF in read_geos5fcst')
     call LIS_verify(nf90_get_var(ftn,swlandid,swland),&
          'nf90_get_var failed for SWLAND in read_geos5fcst')
!EB

     call LIS_verify(nf90_close(ftn),&
          'nf90_open failed in read_geos5fcst')

     call interp_geos5fcst(n, findex,tlml, varfield)
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 geos5fcst_struc(n)%metdata1(1,m,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 geos5fcst_struc(n)%metdata2(1,m,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo

     call interp_geos5fcst(n, findex,qlml, varfield)
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 geos5fcst_struc(n)%metdata1(2,m,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 geos5fcst_struc(n)%metdata2(2,m,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo

     call interp_geos5fcst(n, findex,swgdn, varfield)
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 geos5fcst_struc(n)%metdata1(3,m,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 geos5fcst_struc(n)%metdata2(3,m,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo

     call interp_geos5fcst(n, findex,lwgab, varfield)
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 geos5fcst_struc(n)%metdata1(4,m,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 geos5fcst_struc(n)%metdata2(4,m,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo

     call interp_geos5fcst(n, findex,speedlml, varfield)

     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 geos5fcst_struc(n)%metdata1(5,m,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
                 geos5fcst_struc(n)%metdata1(6,m,&
                      LIS_domain(n)%gindex(c,r)) &
                      = 0.0
              elseif(order.eq.2) then 
                 geos5fcst_struc(n)%metdata2(5,m,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
                 geos5fcst_struc(n)%metdata2(6,m,&
                      LIS_domain(n)%gindex(c,r)) &
                      = 0.0
              endif
           endif
        end do
     enddo


     call interp_geos5fcst(n, findex,ps, varfield)
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 geos5fcst_struc(n)%metdata1(7,m,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 geos5fcst_struc(n)%metdata2(7,m,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo

     do r=1,geos5fcst_struc(n)%nr
        do c=1,geos5fcst_struc(n)%nc
           prectot(c,r)= prectot(c,r) - precsno(c,r)
           if(prectot(c,r).lt.0) then 
              prectot(c,r) = precsno(c,r)
           endif
        enddo
     enddo

     call interp_geos5fcst(n, findex,prectot, varfield)
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 geos5fcst_struc(n)%metdata1(8,m,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 geos5fcst_struc(n)%metdata2(8,m,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo

     call interp_geos5fcst(n, findex,precsno, varfield)
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 geos5fcst_struc(n)%metdata1(9,m,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 geos5fcst_struc(n)%metdata2(9,m,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo

     call interp_geos5fcst(n, findex,preccon, varfield)
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 geos5fcst_struc(n)%metdata1(10,m,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 geos5fcst_struc(n)%metdata2(10,m,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo
!SB
     call interp_geos5fcst(n, findex,hlml, varfield)
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then
              if(order.eq.1) then
                 geos5fcst_struc(n)%metdata1(11,m,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then
                 geos5fcst_struc(n)%metdata2(11,m,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo

     call interp_geos5fcst(n, findex,pardr, varfield)
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then
              if(order.eq.1) then
                 geos5fcst_struc(n)%metdata1(12,m,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then
                 geos5fcst_struc(n)%metdata2(12,m,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo

     call interp_geos5fcst(n, findex,pardf, varfield)
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then
              if(order.eq.1) then
                 geos5fcst_struc(n)%metdata1(13,m,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then
                 geos5fcst_struc(n)%metdata2(13,m,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo


     call interp_geos5fcst(n, findex,swland, varfield)
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then
              if(order.eq.1) then
                 geos5fcst_struc(n)%metdata1(14,m,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then
                 geos5fcst_struc(n)%metdata2(14,m,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
           endif
        end do
     enddo


!EB

  else
     write(LIS_logunit,*) &
          '[ERR] Could not find file: ',trim(name)
     ferror = 0
  endif
#endif     

end subroutine read_geos5fcst


