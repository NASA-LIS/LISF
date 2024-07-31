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
! !ROUTINE: read_HiMATGMU
! \label{read_HiMATGMU}
!
! !REVISION HISTORY:
!  28 Jul 2017: Sujay Kumar;  Data and code implementation
!
! !INTERFACE:
subroutine read_HiMATGMU( n, fname,findex,order, ferror_HiMATGMU )

! !USES:
  use LIS_coreMod,        only : LIS_rc, LIS_domain
  use LIS_logMod,         only : LIS_logunit, LIS_verify
  use LIS_metforcingMod,  only : LIS_forc
  use HiMATGMU_forcingMod,    only : HiMATGMU_struc

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  character(len=*)   :: fname          
  integer, intent(in) :: findex
  integer, intent(in) :: order
  integer             :: ferror_HiMATGMU

! !DESCRIPTION:
!  For the given time, reads parameters from STAGE4 datasets
!  and interpolates to a designated user-domain.
!  NOTE:: These subroutines use the READ\_GRIB routines for
!         for opening and reading the STAGE IV grib files.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[fname]
!    name of the hourly STAGE4 file
!  \item[ferror\_HiMATGMU]
!    flag to indicate success of the call (=0 indicates success)
!  \item[filehr]
!    current file hour
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_HiMATGMU](\ref{interp_HiMATGMU}) \newline
!    spatially interpolates the STAGE4 data
!  \end{description}
!EOP

  integer     :: i, j, t, iret ! Loop indicies and error flags
  integer     :: ftn, ndata, index1
  real        :: precip_regrid(LIS_rc%lnc(n),LIS_rc%lnr(n))     
  real        :: pcpin(HiMATGMU_struc(n)%ncol, HiMATGMU_struc(n)%nrow)
  real        :: HiMATGMUin(HiMATGMU_struc(n)%ncol*HiMATGMU_struc(n)%nrow)
  logical*1   :: lb(HiMATGMU_struc(n)%ncol*HiMATGMU_struc(n)%nrow)
  logical     :: file_exists            ! Check Stage IV file status 
  integer     :: pcpId,c,r
  real        :: missingValue

!=== End Variable Definition =======================
  if(order.eq.1) then 
     HiMATGMU_struc(n)%metdata1 = LIS_rc%udef
  elseif(order.eq.2) then 
     HiMATGMU_struc(n)%metdata2 = LIS_rc%udef
  endif

!-- Set necessary parameters for call to interp_HiMATGMU   
  precip_regrid = LIS_rc%udef

!-- Check initially if file exists:
  inquire (file=fname, exist=file_exists ) ! Check if file exists
  if (.not. file_exists)  then 
     write(LIS_logunit,*)"** Missing HiMAT GMU precipitation file: ", trim(fname)
     ferror_HiMATGMU = 1
     return
  endif

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

  ndata = HiMATGMU_struc(n)%ncol * HiMATGMU_struc(n)%nrow
  HiMATGMUin = 0.0

  call LIS_verify(nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=ftn),&
       'nf90_open failed in read_HiMATGMU')

  call LIS_verify(nf90_inq_varid(ftn,'PRECTOT_ds',pcpId),&
       'nf90_inq_varid failed for PRECTOT_ds in read_HiMATGMU')

  call LIS_verify(nf90_get_var(ftn,pcpid,pcpin),&
       'nf90_get_var failed for PRECTOT_ds in read_HiMATGMU')     

  call LIS_verify(nf90_close(ftn))
  
  lb = .false.
  do r=1,HiMATGMU_struc(n)%nrow
     do c=1,HiMATGMU_struc(n)%ncol
        HiMATGMUin(c+(r-1)*HiMATGMU_struc(n)%ncol) = & 
             pcpin(c,r)
        if (pcpin(c,r) .ne. missingvalue ) then
           lb(c+(r-1)*HiMATGMU_struc(n)%ncol) = .true. 
        endif
     enddo
  enddo

  call interp_HiMATGMU( n, findex, ndata, &
       HiMATGMUin, lb, LIS_rc%gridDesc(n,:), &
       LIS_rc%lnc(n), LIS_rc%lnr(n), precip_regrid )
  
  do j = 1, LIS_rc%lnr(n)
     do i = 1, LIS_rc%lnc(n)
        
        !if ( precip_regrid(i,j) .ne. -1.0 ) then
        if ( precip_regrid(i,j) .ne. LIS_rc%udef ) then
           index1 = LIS_domain(n)%gindex(i,j)
           if(index1 .ne. -1) then
              if(order.eq.1) then 
                 HiMATGMU_struc(n)%metdata1(1,index1) = precip_regrid(i,j) 
              elseif(order.eq.2) then 
                 HiMATGMU_struc(n)%metdata2(1,index1) = precip_regrid(i,j) 
              endif
           endif
        endif
        
     enddo
  enddo

#endif
     
end subroutine read_HiMATGMU

