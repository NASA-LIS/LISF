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
! !ROUTINE: read_gdasbc
!  \label{read_gdasbc}
!
! !REVISION HISTORY:
! 
! !INTERFACE:
subroutine read_gdasbc(n, kk, findex, order, month, name,ferror)
! !USES:
  use LIS_coreMod
  use LIS_logMod, only         : LIS_logunit, LIS_verify, LIS_warning
  use LIS_metforcingMod, only  : LIS_forc
  use gdasbc_forcingMod, only  : gdasbc_struc
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
!  GDASBC data, transforms into 11 LIS forcing 
!  parameters and interpolates to the LIS domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous 
!    hourly instance, order=2, read the next hourly instance)
!  \item[n]
!    index of the nest
!  \item[name]
!    name of the hourly GDASBC forecast file
!  \item[ferror]
!    flag to indicate success of the call (=1 indicates success)
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_gdasbc](\ref{interp_gdasbc}) \newline
!    spatially interpolates a GDASBC variable
!  \end{description}
!EOP

  integer                   :: iv, ftn
  integer                   :: gdasbc
  integer                   :: k,t,c,r,iret,rc
  integer                   :: pcpId
  real                      :: missingValue 
  integer                   :: nvars
  integer                   :: igrib
  logical                   :: pcp_flag, var_found
  logical                   :: var_status(11)
  logical                   :: file_exists
  logical*1                 :: lb(gdasbc_struc(n)%ncold*gdasbc_struc(n)%nrold)
  real                      :: pcp(gdasbc_struc(n)%ncold,gdasbc_struc(n)%nrold)
  real                      :: pcpin(gdasbc_struc(n)%ncold*gdasbc_struc(n)%nrold)
  real                      :: varfield(LIS_rc%lnc(n),LIS_rc%lnr(n))

  ferror = 1
  iv = 0
#if (defined USE_NETCDF3 || defined USE_NETCDF4)

  inquire (file=name, exist=file_exists)
  if (file_exists) then      
  
     call LIS_verify(nf90_open(path=trim(name), mode=NF90_NOWRITE, &
          ncid=ftn), 'nf90_open failed for in read_gdasbc')
     
     call LIS_verify(nf90_inq_varid(ftn,'TotalPrecip_tavg',pcpId), &
          'nf90_inq_varid failed for TotalPrecip_tavg in read_gdasbc')
     
     call LIS_verify(nf90_get_var(ftn,pcpId, pcp), &
          'nf90_get_var failed for TotalPrecip_tavg in read_gdasbc')
     call LIS_verify(nf90_close(ftn))
     
     pcp_flag = .true.

     lb = .false. 
     do r=1,gdasbc_struc(n)%nrold
        do c=1,gdasbc_struc(n)%ncold
           pcpin(c+(r-1)*gdasbc_struc(n)%ncold) = &
                pcp(c,r)
           if(pcp(c,r).ge.0) then
              lb(c+(r-1)*gdasbc_struc(n)%ncold) = .true.
           endif
        enddo
     enddo

           
     call interp_gdasbc(n, findex, month, pcp_flag, gdasbc,&
          pcpin,&
          lb, LIS_rc%gridDesc(n,:), &
          LIS_rc%lnc(n),LIS_rc%lnr(n),varfield )

     
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 gdasbc_struc(n)%metdata1(kk,1,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 gdasbc_struc(n)%metdata2(kk,1,&
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

end subroutine read_gdasbc


!BOP
! !ROUTINE: interp_gdasbc
! \label{interp_gdasbc}
!
! !INTERFACE:
subroutine interp_gdasbc(n,findex, month, pcp_flag, &
     input_size,input_data,input_bitmap,&
     lis_gds,nc,nr, &
     output_2d)
! !USES:
  use LIS_coreMod, only : LIS_rc, LIS_domain
  use gdasbc_forcingMod, only :gdasbc_struc
  
  implicit none

! !ARGUMENTS:   
  integer, intent(in)   :: n 
  integer, intent(in)   :: findex
  integer, intent(in)   :: month
  logical, intent(in)   :: pcp_flag
  integer, intent(in)   :: input_size
  real, intent(in)      :: input_data(input_size)
  logical*1, intent(in) :: input_bitmap(input_size)
  real, intent(in)      :: lis_gds(50)
  integer, intent(in)   :: nc
  integer, intent(in)   :: nr
  real, intent(inout)   :: output_2d(nc,nr)
!
! !DESCRIPTION:
!   This subroutine interpolates a given GDASBC field 
!   to the LIS grid. 
!  The arguments are: 
!  \begin{description}
! \item[n]
!  index of the nest
! \item[kpds]
!  grib decoding array
! \item[input\_size]
!  number of elements in the input grid
! \item[f]
!  input data array to be interpolated
! \item[input\_bitmap]
!  input bitmap
! \item[lis\_gds]
!  array description of the LIS grid
! \item[nc]
!  number of columns (in the east-west dimension) in the LIS grid
! \item[nr]
!  number of rows (in the north-south dimension) in the LIS grid
! \item[output\_2d]
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
!  \item[neighbor\_interp](\ref{neighbor_interp}) \newline
!    spatially interpolate the forcing data using nearest neighbor interpolation
! \end{description}
!EOP
  integer :: iret
  integer :: mo
  integer :: count1,i,j
  real, dimension(nc*nr) :: output_data
  logical*1 :: output_bitmap(nc*nr)

!=== End variable declarations

  mo = nc*nr
!-----------------------------------------------------------------------
! Initialize output bitmap. 
!-----------------------------------------------------------------------
  output_bitmap = .true.
!-----------------------------------------------------------------------  
! Interpolate to LIS grid
!-----------------------------------------------------------------------
  select case( LIS_rc%met_interp(findex) )

    case( "bilinear" )
     call bilinear_interp(lis_gds,input_bitmap,input_data,&
          output_bitmap,&
          output_data,gdasbc_struc(n)%mi,mo,&
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          gdasbc_struc(n)%w111, gdasbc_struc(n)%w121,&
          gdasbc_struc(n)%w211,gdasbc_struc(n)%w221,&
          gdasbc_struc(n)%n111,gdasbc_struc(n)%n121,&
          gdasbc_struc(n)%n211,gdasbc_struc(n)%n221,LIS_rc%udef,iret)

    case( "budget-bilinear" )
     if (pcp_flag) then     
        call conserv_interp(lis_gds,input_bitmap,input_data,&
             output_bitmap,&
             output_data,gdasbc_struc(n)%mi,mo,& 
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             gdasbc_struc(n)%w112,gdasbc_struc(n)%w122,&
             gdasbc_struc(n)%w212,gdasbc_struc(n)%w222,&
             gdasbc_struc(n)%n112,gdasbc_struc(n)%n122,&
             gdasbc_struc(n)%n212,gdasbc_struc(n)%n222,LIS_rc%udef,iret)
     else 
        call bilinear_interp(lis_gds,input_bitmap,input_data,&
             output_bitmap,&
             output_data,gdasbc_struc(n)%mi,mo,&
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             gdasbc_struc(n)%w111,gdasbc_struc(n)%w121,&
             gdasbc_struc(n)%w211,gdasbc_struc(n)%w221,&
             gdasbc_struc(n)%n111,gdasbc_struc(n)%n121,&
             gdasbc_struc(n)%n211,gdasbc_struc(n)%n221,LIS_rc%udef,iret)
     endif

    case( "neighbor" )
     call neighbor_interp(lis_gds,input_bitmap,input_data,&
          output_bitmap,&
          output_data,gdasbc_struc(n)%mi,mo,&
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          gdasbc_struc(n)%n113,LIS_rc%udef,iret)

  end select


!-----------------------------------------------------------------------    
! convert the interpolated data to 2d. 
!-----------------------------------------------------------------------    
  count1 = 0
  do j = 1, nr
     do i = 1, nc
        output_2d(i,j) = output_data(i+count1)
     enddo
     count1 = count1 + nc
  enddo

end subroutine interp_gdasbc
