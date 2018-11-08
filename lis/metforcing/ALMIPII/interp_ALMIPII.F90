!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: interp_ALMIPII
! \label{interp_ALMIPII}
!
! !INTERFACE:
subroutine interp_ALMIPII(n, findex,fvar,varfield,precip_field)
! !USES:
  use LIS_coreMod
  use ALMIPII_forcingMod,  only : ALMIPII_struc
  
  implicit none

! !ARGUMENTS:   
  integer, intent(in)   :: n 
  integer, intent(in)   :: findex
  real, intent(in)      :: fvar(ALMIPII_struc(n)%nc,ALMIPII_struc(n)%nr)
  real, intent(inout)   :: varfield(LIS_rc%ngrid(n))
  logical               :: precip_field
!
! !DESCRIPTION:
!   This subroutine interpolates a given ALMIPII field 
!   to the LIS grid. 
!  The arguments are: 
!  \begin{description}
! \item[n]
!  index of the nest
! \item[f]
!  input data array to be interpolated
! \item[lb]
!  input bitmap
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
!  \item[neighbor\_interp](\ref{neighbor_interp}) \newline
!    spatially interpolate the forcing data using neighbor interpolation
! \end{description}
!EOP
!  real    :: f(ALMIPII_struc(n)%mi)
!  logical*1 :: lb(ALMIPII_struc(n)%mi)
!  integer :: ip, ibi,km,iret
!  integer :: ibo,mo
  real    :: array(LIS_rc%lnc(n),LIS_rc%lnr(n))
  integer :: count1,i,j
  real    :: fvar_t(ALMIPII_struc(n)%nc,ALMIPII_struc(n)%nr)
!  real, dimension(LIS_rc%lnc(n)*LIS_rc%lnr(n)) :: lis1d
!  logical*1 :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))

  do j=1,LIS_rc%gnr(n)
     do i=1,LIS_rc%gnc(n)
!        fvar_t(i,j) = fvar(i,LIS_rc%gnr(n)-j+1)
        fvar_t(i,j) = fvar(i,j)
     enddo
  enddo

  array(:,:) = &
       fvar_t(LIS_ews_halo_ind(n,LIS_localPet+1):&         
       LIS_ewe_halo_ind(n,LIS_localPet+1), &
       LIS_nss_halo_ind(n,LIS_localPet+1): &
       LIS_nse_halo_ind(n,LIS_localPet+1))

  do j = 1, LIS_rc%lnr(n)
     do i = 1, LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(i,j).ne.-1) then
           varfield(LIS_domain(n)%gindex(i,j)) = array(i,j)
        endif
     enddo
  enddo

#if 0 
!=== End variable declarations

  lb = .true. 

  do i=1,ALMIPII_struc(n)%nc
     do j=1,ALMIPII_struc(n)%nr
        f(i+(j-1)*ALMIPII_struc(n)%nc) = fvar(i,j)
        if(fvar(i,j).gt.1E10) lb(i+(j-1)*ALMIPII_struc(n)%nc) = .false.
     enddo
  enddo

!-----------------------------------------------------------------------
! Setting interpolation options (ip=0,bilinear)
! (km=1, one parameter, ibi=1,use undefined bitmap
! (needed for soil moisture and temperature only)
! Use budget bilinear (ip=3) for precip forcing fields
!-----------------------------------------------------------------------
  mo = LIS_rc%lnc(n)*LIS_rc%lnr(n)
  km=1

!-----------------------------------------------------------------------
! Initialize output bitmap. Important for soil moisture and temp.
!-----------------------------------------------------------------------
  lo = .true.

!-----------------------------------------------------------------------  
! Interpolate to LIS grid
!-----------------------------------------------------------------------
  if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 
     call bilinear_interp(LIS_rc%gridDesc(n,:),lb,f,lo,&
          lis1d,ALMIPII_struc(n)%mi,mo,&
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          ALMIPII_struc(n)%w111, ALMIPII_struc(n)%w121,&
          ALMIPII_struc(n)%w211,ALMIPII_struc(n)%w221,&
          ALMIPII_struc(n)%n111,ALMIPII_struc(n)%n121,&
          ALMIPII_struc(n)%n211,ALMIPII_struc(n)%n221,LIS_rc%udef,iret)
  elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 
     if(precip_field) then 
        call conserv_interp(LIS_rc%gridDesc(n,:),lb,f,lo,&
             lis1d,ALMIPII_struc(n)%mi,mo,& 
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             ALMIPII_struc(n)%w112,ALMIPII_struc(n)%w122,&
             ALMIPII_struc(n)%w212,ALMIPII_struc(n)%w222,&
             ALMIPII_struc(n)%n112,ALMIPII_struc(n)%n122,&
             ALMIPII_struc(n)%n212,ALMIPII_struc(n)%n222,LIS_rc%udef,iret)
     else 
        call bilinear_interp(LIS_rc%gridDesc(n,:),lb,f,lo,&
             lis1d,ALMIPII_struc(n)%mi,mo,&
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             ALMIPII_struc(n)%w111,ALMIPII_struc(n)%w121,&
             ALMIPII_struc(n)%w211,ALMIPII_struc(n)%w221,&
             ALMIPII_struc(n)%n111,ALMIPII_struc(n)%n121,&
             ALMIPII_struc(n)%n211,ALMIPII_struc(n)%n221,LIS_rc%udef,iret)
     endif
  elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then 
        call neighbor_interp(LIS_rc%gridDesc(n,:),lb,f,lo,&
             lis1d,ALMIPII_struc(n)%mi,mo,&
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             ALMIPII_struc(n)%n113,LIS_rc%udef,iret)
  endif
!-----------------------------------------------------------------------    
! Create 2D array for main program. Also define a "soil" mask
! due to different geography between LDAS & LDAS. For LDAS land 
! points not included in LDAS geography dataset only.
!-----------------------------------------------------------------------    
  count1 = 0

!     open(100,file='almip_ip.bin',form='unformatted')
!     write(100) lis1d
!     close(100)
!     stop
#endif
end subroutine interp_ALMIPII
