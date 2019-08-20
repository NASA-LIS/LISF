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
! !ROUTINE: read_capa
! \label{read_capa}
!
! !REVISION HISTORY:
! 23 Nov 2009: Sujay Kumar, Initial Specification
! 
! !INTERFACE:
subroutine read_capa(n, findex, fname, ferror_capa)
! !USES:
  use LIS_coreMod,        only : LIS_rc, LIS_domain
  use LIS_logMod,         only : LIS_logunit, LIS_endrun
  use LIS_metforcingMod,  only : LIS_forc
  use capa_forcingMod,    only : capa_struc

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex
  character(len=80)   :: fname          
  integer             :: ferror_capa
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  CAPA data and interpolates to the LIS domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!   \item[findex]
!    index of the forcing source
!  \item[fname]
!    name of the 6 hour CAPA file
!  \item[ferror\_capa]
!    flag to indicate success of the call (=0 indicates success)
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_capa](\ref{interp_capa}) \newline
!    spatially interpolates the CAPA data
!  \end{description}
!EOP

#if (defined USE_GRIBAPI)
  integer               :: i,j,iret,jret  ! Loop indicies and error flags
  real, allocatable         :: precip_regrid(:,:)  ! Interpolated precipitation array
  integer               :: ncapa
  real, allocatable     :: capain(:)
  logical*1,allocatable :: lb(:)
  integer               :: ftn, igrib, index1
!=== End Variable Definition =======================

  allocate (precip_regrid(LIS_rc%lnc(n),LIS_rc%lnr(n)))

  precip_regrid = -1.0    
!------------------------------------------------------------------------    
! Set necessary parameters for call to interp_gdas    
!------------------------------------------------------------------------    
  ncapa = capa_struc(n)%ncold*capa_struc(n)%nrold
  allocate(capain(ncapa))
  allocate(lb(ncapa)) 
  capain = 0.0

  call grib_open_file(ftn,trim(fname),'r',iret)
  if (iret == 0 ) then
     ferror_capa = 1 ! assume successful read

     ! There is only one message per input file.
     call grib_new_from_file(ftn, igrib, iret)
     call grib_get(igrib,'values',capain,iret)
     if(iret /= 0) then 
        ferror_capa = 0
        write(LIS_logunit,*) 'ERR: read_capa: ' //     &
                             'Could not retrieve values in file: ' // &
                             trim(fname)
     endif
     call grib_release(igrib,iret)
     call grib_close_file(ftn)

     lb = .true.
     do i=1,ncapa
        if(capain(i).ge.52.648) then 
           lb(i) = .false. 
        endif
     enddo

     call interp_capa(n,ncapa,capain,lb,LIS_rc%gridDesc, &
          LIS_rc%lnc(n),LIS_rc%lnr(n),precip_regrid)

     do j = 1,LIS_rc%lnr(n)
        do i = 1,LIS_rc%lnc(n)
           if (precip_regrid(i,j) .ne. -1.0) then
              index1 = LIS_domain(n)%gindex(i,j)
              if(index1 .ne. -1) then 
                 capa_struc(n)%metdata2(1,index1) = precip_regrid(i,j)
              endif
           endif
        enddo
     enddo
     
     write(LIS_logunit,*) "Obtained CAPA CPC precipitation data ", trim(fname)
  else
     ferror_capa = 0
     write(LIS_logunit,*) "Missing CAPA CPC precipitation data ", trim(fname)
  endif
  deallocate (precip_regrid)
  deallocate(lb)
  deallocate(capain)
#else
  write(LIS_logunit,*) 'ERR: read_capa: ' //              &
                       'CaPA support requires the GRIB_API library. ', &
                       'Please recompile LIS.'

  call LIS_endrun
#endif
end subroutine read_capa

