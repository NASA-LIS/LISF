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
! !ROUTINE: read_nldas2a
!  \label{read_nldas2a}
!
! !REVISION HISTORY:
!  11 Apr 2000: Brian Cosgrove; changed code to use Forcing Mask (With
!               inland water filled in).  Deleted unused variables.
!  27 Apr 2000: Brian Cosgrove; changed code to use the original 
!               mask again since that is the
!               mask which  NCEP has already applied to the forcing data
!               by the time NASA gets it......not possible to use the 
!               expanded NASA forcing mask
!  1  May 2000: Brian Cosgrove; changed code so that if parameter 11 (sw)
!               is not found in hourly ncep data, it will just use
!               edas-based shortwave from the hourly ncep files
!  20 Jun 2000: Brian Cosgrove; changed code so that it uses  LDAS%UDEF and
!               not a hard-wired undefined value of -999.9 and -999.0
!  18 Aug 2000: Brian Cosgrove; changed code so that FMASK and not MASK
!               is used when ungribbing.  NCEP data already has a mask applied
!               to it and so may not be able to supply forcing data to
!               all LDAS land forcing points.  In areas where LDAS
!               forcing mask states that land exists, but where NCEP forcing
!               data is non-existant, assign undefined value to forcing data.
!  22 Aug 2000: Brian Cosgrove; Altered code for US/Mexico/Canada Mask
!  05 Sep 2001: Brian Cosgrove; Removed dirnom and infile variables, changed
!               call to ungribncep to match removal.  Added code to make use
!               of precip weighting mask
!  02 Feb 2004: Sujay Kumar; Initial Specification in LDT
!  24 Aug 2007: Chuck Alonge; Modified for use with NLDAS2 data
!  25 Jan 2012: Sujay Kumar; Switched to the use of grib-api library
!  14 Mar 2014: David Mocko: Added CAPE and PET forcing from NLDAS-2
! 
! !INTERFACE:
subroutine read_nldas2a(n, findex, order, month, name,ferror)
! !USES:
  use LDT_coreMod
  use LDT_logMod, only         : LDT_logunit, LDT_verify, LDT_warning
  use LDT_metforcingMod, only  : LDT_forc
  use nldas2_forcingMod, only  : nldas2_struc
#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in)      :: findex
  integer, intent(in)      :: n
  integer, intent(in)      :: order
  integer, intent(out)     :: month
  character(len=*), intent(in) :: name
  integer, intent(out)     :: ferror
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  NLDAS2 data, transforms into 11 LDT forcing 
!  parameters and interpolates to the LDT domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous 
!    hourly instance, order=2, read the next hourly instance)
!  \item[n]
!    index of the nest
!  \item[name]
!    name of the hourly NLDAS2 forecast file
!  \item[ferror]
!    flag to indicate success of the call (=1 indicates success)
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_nldas2](\ref{interp_nldas2}) \newline
!    spatially interpolates a NLDAS2 variable
!  \end{description}
!EOP

  integer                   :: iv, ftn
  integer                   :: nldas2
  integer                   :: kk,t,c,r,iret,rc
  integer                   :: var_index
  real                      :: missingValue 
  integer                   :: nvars
  integer                   :: igrib
  logical                   :: pcp_flag, var_found
  logical                   :: var_status(11)
  logical                   :: file_exists
  integer                   :: pds5(11), pds7(11), pds2(11)
  integer                   :: pds5_val, pds7_val, pds2_val
  logical*1, allocatable    :: lb(:)
  real, allocatable         :: f(:)
  real, allocatable         :: nldas_forcing(:,:)
  real                      :: varfield(LDT_rc%lnc(n),LDT_rc%lnr(n))

  ferror = 1
  iv = 0
#if(defined USE_GRIBAPI) 
  nldas2 = (nldas2_struc(n)%nc*nldas2_struc(n)%nr)
! Order of variables to be assigned in the NLDAS-2 Forcing "A" files.
! Note that this is NOT the order of the fields in the actual files,
! but the order in which they are assigned to "metdata". - dmm
  pds5 = (/ 011,051,204,205,033,034,001,061,153,228,157/) !parameter
  pds7 = (/ 002,002,000,000,010,010,000,000,000,000,180/) !htlev2
  pds2 = (/  84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84/)
  
  allocate(nldas_forcing(nldas2_struc(n)%nc*nldas2_struc(n)%nr,11))

  varfield = 0 
  ferror = 1
  var_status = .false. 

  inquire (file=name, exist=file_exists)
  if (file_exists) then      

     call grib_open_file(ftn,trim(name),'r',iret)
     if(iret.ne.0) then 
        write(LDT_logunit,*) &
             'Could not open file: ',trim(name)
        ferror = 0
        return
     endif

     call grib_count_in_file(ftn,nvars,iret)
     call LDT_verify(iret, 'error in grib_count_in_file in read_nldas2a')

     allocate(lb(nldas2_struc(n)%nc*nldas2_struc(n)%nr))
     allocate(f(nldas2_struc(n)%nc*nldas2_struc(n)%nr))
     
     do kk=1,nvars
        call grib_new_from_file(ftn, igrib, iret)
        call LDT_warning(iret, 'error in grib_new_from_file in read_nldas2a')
        if(iret.ne.0) then 
           write(LDT_logunit,*) &
                'Could not retrieve entries in file: ',trim(name)
           ferror = 0
           deallocate(lb)
           deallocate(f)
           return           
        endif

        call grib_get(igrib,'indicatorOfParameter',pds5_val,rc)
        call LDT_verify(rc, 'error in grib_get: indicatorOfParameter in read_nldas2a')

        call grib_get(igrib,'level',pds7_val,rc)
        call LDT_verify(rc, 'error in grib_get: level in read_nldas2a')

        call grib_get(igrib,'generatingProcessIdentifier',pds2_val,rc)
        call LDT_verify(rc, 'error in grib_get: generatingProcessIdentifier in read_nldas2a')

        var_found = .false. 
        do iv=1,nvars
           if((pds5_val.eq.pds5(iv)).and.&
                (pds7_val.eq.pds7(iv)).and.&
                (pds2_val.eq.pds2(iv))) then
              var_found = .true.
              var_index = iv
              var_status(iv) = .true. 
              exit
           endif
        enddo
        f = -9999.0
        call grib_get(igrib,'values',f,rc)
        call LDT_warning(rc, 'error in grib_get:values in read_nldas2a')

        if(rc.ne.0) then 
           write(LDT_logunit,*) &
                'Could not retrieve entries in file: ',trim(name)
           ferror = 0
           deallocate(lb)
           deallocate(f)
           return           
        endif

        call grib_get(igrib,'missingValue',missingValue,rc)
        call LDT_verify(rc, 'error in grib_get:missingValue in read_nldas2a')

        call grib_release(igrib,rc)
        call LDT_verify(rc, 'error in grib_release in read_nldas2a')
        
        if(var_found) then 
           lb = .false.
           do t=1,nldas2
              if(f(t).ne.missingValue) then 
                 nldas_forcing(t,var_index) = f(t)
                 lb(t) = .true. 
              else
                 nldas_forcing(t,var_index) = LDT_rc%udef
              endif
           enddo
        endif
     enddo

     call grib_close_file(ftn)
     deallocate(f)     
     
     do kk=1,nvars
        if(.not.var_status(kk)) then 
           write(LDT_logunit,*) &
                'Could not retrieve entries in file: ',trim(name)
           ferror = 0
           
           return
        endif
     enddo

     if(trim(LDT_rc%met_ecor(findex)).ne."none") then 
        do t=1,nldas2
           if(lb(t)) then 
              call nldas2_ec_removal(n,t, &
                   nldas_forcing(t,1), &
                   nldas_forcing(t,2), &
                   nldas_forcing(t,4), &
                   nldas_forcing(t,7))
           endif
        enddo
     endif

     do iv = 1,nvars
        pcp_flag = .false. 
        if(iv.eq.8.or.iv.eq.9) pcp_flag = .true.            
        
        call interp_nldas2(n, findex,month,pcp_flag,nldas2,&
             nldas_forcing(:,iv), & 
             lb,LDT_rc%gridDesc(n,:), &
             LDT_rc%lnc(n),LDT_rc%lnr(n),varfield)
        
        do r=1,LDT_rc%lnr(n)
           do c=1,LDT_rc%lnc(n)
              if(LDT_domain(n)%gindex(c,r).ne.-1) then 
                 if(order.eq.1) then 
                    LDT_forc(n,findex)%metdata1(iv,&
                         LDT_domain(n)%gindex(c,r)) &
                         = varfield(c,r)
                 elseif(order.eq.2) then 
                    LDT_forc(n,findex)%metdata2(iv,&
                         LDT_domain(n)%gindex(c,r))&
                         = varfield(c,r)
                 endif
              endif
           end do
        enddo
     enddo

     deallocate(lb)
     deallocate(nldas_forcing)    
  else
     write(LDT_logunit,*) &
          'Could not find file: ',trim(name)
     ferror = 0
  endif
#endif     

end subroutine read_nldas2a


!BOP
! !ROUTINE: interp_nldas2
! \label{interp_nldas2}
!
! !INTERFACE:
subroutine interp_nldas2(n,findex, month, pcp_flag, &
     input_size,input_data,input_bitmap,&
     ldt_gds, nc, nr, output_2d)

! !USES:
  use LDT_coreMod, only : LDT_rc, LDT_domain
  use LDT_spatialDownscalingMod
  use nldas2_forcingMod, only :nldas2_struc
  
  implicit none

! !ARGUMENTS:   
  integer, intent(in)   :: n 
  integer, intent(in)   :: findex
  integer, intent(in)   :: month
  logical, intent(in)   :: pcp_flag
  integer, intent(in)   :: input_size
  real, intent(in)      :: input_data(input_size)
  logical*1, intent(in) :: input_bitmap(input_size)
  real, intent(in)      :: ldt_gds(20)
  integer, intent(in)   :: nc
  integer, intent(in)   :: nr
  real, intent(inout)   :: output_2d(nc,nr)
!
! !DESCRIPTION:
!   This subroutine interpolates a given NLDAS2 field 
!   to the LDT grid. 
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
! \item[ldt\_gds]
!  array description of the LDT grid
! \item[nc]
!  number of columns (in the east-west dimension) in the LDT grid
! \item[nr]
!  number of rows (in the north-south dimension) in the LDT grid
! \item[output\_2d]
!  output interpolated field
!  \end{description} 
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
! Apply downscaling
!-----------------------------------------------------------------------    

#if 0
  if(pcp_flag.and.LDT_rc%pcp_downscale(findex).ne.0) then 
!input_data becomes the ratio field. 
     call LDT_generatePcpClimoRatioField(n,findex,"NLDAS2",&
          month, & 
          input_size, &
          input_data, &
          input_bitmap)     
  endif
#endif

!-----------------------------------------------------------------------  
! Interpolate to LDT grid
!-----------------------------------------------------------------------
  select case( LDT_rc%met_gridtransform(findex) )
   case( "bilinear" )
     call bilinear_interp(ldt_gds,input_bitmap,input_data,&
          output_bitmap,&
          output_data,nldas2_struc(n)%mi,mo,&
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          nldas2_struc(n)%w111, nldas2_struc(n)%w121,&
          nldas2_struc(n)%w211,nldas2_struc(n)%w221,&
          nldas2_struc(n)%n111,nldas2_struc(n)%n121,&
          nldas2_struc(n)%n211,nldas2_struc(n)%n221,LDT_rc%udef,iret)
   case( "budget-bilinear" )  
     if (pcp_flag) then     
        call conserv_interp(ldt_gds,input_bitmap,input_data,&
             output_bitmap,&
             output_data,nldas2_struc(n)%mi,mo,& 
             LDT_domain(n)%lat, LDT_domain(n)%lon,&
             nldas2_struc(n)%w112,nldas2_struc(n)%w122,&
             nldas2_struc(n)%w212,nldas2_struc(n)%w222,&
             nldas2_struc(n)%n112,nldas2_struc(n)%n122,&
             nldas2_struc(n)%n212,nldas2_struc(n)%n222,LDT_rc%udef,iret)
     else 
        call bilinear_interp(ldt_gds,input_bitmap,input_data,&
             output_bitmap,&
             output_data,nldas2_struc(n)%mi,mo,&
             LDT_domain(n)%lat, LDT_domain(n)%lon,&
             nldas2_struc(n)%w111,nldas2_struc(n)%w121,&
             nldas2_struc(n)%w211,nldas2_struc(n)%w221,&
             nldas2_struc(n)%n111,nldas2_struc(n)%n121,&
             nldas2_struc(n)%n211,nldas2_struc(n)%n221,LDT_rc%udef,iret)
     endif
   case( "neighbor" )
     call neighbor_interp(ldt_gds,input_bitmap,input_data,&
          output_bitmap,&
          output_data,nldas2_struc(n)%mi,mo,&
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          nldas2_struc(n)%n113,LDT_rc%udef,iret)
  end select

#if 0
  if(pcp_flag.and.LDT_rc%pcp_downscale(findex).ne.0) then 

     call LDT_pcpClimoDownscaling(n,findex,month,&
          nc*nr, output_data, output_bitmap)
  endif
#endif

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

end subroutine interp_nldas2
