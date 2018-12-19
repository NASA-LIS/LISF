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
! !ROUTINE: read_AWRAL
! \label{read_AWRAL}
!
! !REVISION HISTORY:
!  30 Jan 2017: Sujay Kumar, Initial version
!
! !INTERFACE:
subroutine read_AWRAL( n, fname,findex,order, ferror_AWRAL )

! !USES:
  use LIS_coreMod,        only : LIS_rc, LIS_domain
  use LIS_logMod,         only : LIS_logunit, LIS_verify
  use LIS_metforcingMod,  only : LIS_forc
  use AWRAL_forcingMod,    only : AWRAL_struc

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  character(len=80)   :: fname          
  integer, intent(in) :: findex
  integer, intent(in) :: order
  integer             :: ferror_AWRAL

! !DESCRIPTION:
!  For the given time, reads parameters from AWRAL datasets
!  and interpolates to a designated user-domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[fname]
!    name of the AWRAL file : file naming convention is currently: AWRALdir/AWRAL_calendaryear_calendarday/var.grb'
!  \item[ferror\_AWRAL]
!    flag to indicate success of the call (=0 indicates success)
!  \item[filehr]
!    current file hour
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_AWRAL](\ref{interp_AWRAL}) \newline
!    spatially interpolates the STAGE4 data
!  \end{description}
!EOP

  integer            :: i, j, t, iret ! Loop indicies and error flags
  integer            :: ftn, ndata, index1,v
  integer, parameter :: N_PF = 6 ! # of AWRAL forcing variables

  character(10) :: var_fname
  real               :: regrid(LIS_rc%lnc(n),LIS_rc%lnr(n))
    character(10), dimension(N_PF), parameter :: awral_fv = (/  &
       'tat       ',    &
       'rgt      ',    &
       'pt     ',    &
       'avpt     ',    &
       'u2t      ',    &
       'radcskyt      '     /)  

  
  real               :: AWRALin(AWRAL_struc(n)%ncol*AWRAL_struc(n)%nrow)
  logical*1          :: lb(AWRAL_struc(n)%ncol*AWRAL_struc(n)%nrow)
  logical            :: file_exists            
  integer            :: igrib
  real               :: missingValue
  integer            :: pds5_val, pds7_val

  integer            :: c1,r1,c,r
  real               :: AWRALin1(AWRAL_struc(n)%ncol*AWRAL_struc(n)%nrow)

!=== End Variable Definition =======================
  if(order.eq.1) then 
     AWRAL_struc(n)%metdata1 = LIS_rc%udef
  elseif(order.eq.2) then 
     AWRAL_struc(n)%metdata2 = LIS_rc%udef
  endif

!-- Set necessary parameters for call to interp_AWRAL   
  !precip_regrid = -1.0
  !pt_regrid = LIS_rc%udef
  !tat_regrid = LIS_rc%udef
  !rgt_regrid = LIS_rc%udef
  !avpt_regrid = LIS_rc%udef
  !u2t_regrid = LIS_rc%udef
  !radcskyt_regrid = LIS_rc%udef

!-- Check initially if file exists:
  do v = 1, N_PF  ! N_AF
    var_fname =awral_fv(v) 
    inquire (file=trim(var_fname), exist=file_exists ) ! Check if file exists
     if (.not. file_exists)  then 
       write(LIS_logunit,*)"[ERR] Missing AWRAL file: ", fname
       ferror_AWRAL = 1
       return
    endif
  enddo

#if(defined USE_GRIBAPI) 
  do v = 1, N_PF ! N_AF
    regrid = LIS_rc%udef
    ndata = AWRAL_struc(n)%ncol * AWRAL_struc(n)%nrow
    AWRALin = 0.0
    var_fname  =awral_fv(v)  
    call grib_open_file(v,trim(var_fname),'r',iret)
    call LIS_verify(iret,'error grib_open_file in read_AWRAL')
  
    call grib_new_from_file(v,igrib,iret)
    call LIS_verify(iret,'error in grib_new_from_file in read_AWRAL')
  
    call grib_get(igrib,'indicatorOfParameter',pds5_val,iret)
    call LIS_verify(iret, 'error in grib_get: indicatorOfParameter in read_AWRAL')
  
    call grib_get(igrib,'level',pds7_val,iret)
    call LIS_verify(iret, 'error in grib_get: level in read_AWRAL')

    if(pds5_val==20 .AND. pds7_val==0) then 
	  call grib_get(igrib,'values',AWRALin,iret)
      call LIS_verify(iret, 'error in grib_get:values in read_AWRAL')
     
    call grib_get(igrib,'missingValue',missingValue,iret)
    call LIS_verify(iret, 'error in grib_get:missingValue in read_AWRAL')
     
    do r=1, AWRAL_struc(n)%nrow
        do c=1,AWRAL_struc(n)%ncol
           c1 = c
           r1 = AWRAL_struc(n)%nrow-r+1
           AWRALin1(c1+(r1-1)*AWRAL_struc(n)%ncol) = &
                AWRALin(c+(r-1)*AWRAL_struc(n)%ncol)
        enddo
    enddo

    lb = .false. 
    do t=1,ndata
        if ( AWRALin1(t) .ne. missingvalue ) then
           lb(t) = .true. 
        endif
    enddo

    call interp_AWRAL( n, findex, ndata, AWRALin1, lb, LIS_rc%gridDesc(n,:), &
          LIS_rc%lnc(n), LIS_rc%lnr(n), regrid )

    do j = 1, LIS_rc%lnr(n)
        do i = 1, LIS_rc%lnc(n)
           
           !if ( precip_regrid(i,j) .ne. -1.0 ) then
           if ( regrid(i,j) .ne. LIS_rc%udef ) then
              index1 = LIS_domain(n)%gindex(i,j)
              if(index1 .ne. -1) then
                 if(order.eq.1) then 
                    AWRAL_struc(n)%metdata1(1,index1) = regrid(i,j) 
                 elseif(order.eq.2) then 
                    AWRAL_struc(n)%metdata2(1,index1) = regrid(i,j) 
                 endif
              endif
           endif
           
        enddo
    enddo

    call grib_release(igrib,iret)
    call LIS_verify(iret,'error in grib_release in read_AWRAL')

  else
    write(LIS_logunit,*) 'Could not retrieve entries in file: ',trim(fname)
    ferror_AWRAL = 1
    return
  endif

  call grib_close_file(v)
  enddo
#endif
     
end subroutine read_AWRAL

