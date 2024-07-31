!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: AGRMET_cmorph
! \label{AGRMET_cmorph}
!
! !REVISION HISTORY:
! 05 May 2013; Moved to AGRMET from suppforcing...Ryan Ruhge/16WS/WXE/SEMS
!
! !INTERFACE:
subroutine AGRMET_cmorph (n, name_cmor, cmorphdata, quad9r )
! !USES:
  use LIS_coreMod, only : LIS_rc, LIS_domain
  use LIS_logMod, only : LIS_logunit, LIS_getNextUnitNumber, &
       LIS_releaseUnitNumber

  implicit none
! !ARGUMENTS:   
  integer, intent(in) :: n
  real                :: cmorphdata(LIS_rc%lnc(n), LIS_rc%lnr(n))
  real, intent(in)    :: quad9r
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  CMORPH data and interpolates to the LIS domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[name\_cmor]
!    name of the CMORPH file
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[AGRMET\_interp\_cmorph](\ref{AGRMET_interp_cmorph}) \newline
!    spatially interpolates the CMORPH data
!  \end{description}
!EOP

  integer :: index1

!==== Local Variables=======================
  integer :: r,ios
  integer :: i,j,xd,yd,ibad
  parameter(xd=4948,yd=1649)    ! Dimension of original CMORPH 8Km data
  integer, parameter :: ncmor=xd*yd
  parameter(ibad=-9999.0)                                  ! Bad (missing data) value
  real  precip(xd,yd)
  real  testout1(xd,yd)                             ! Reconfigured original precip array
  real :: testout(xd,yd)
  real, allocatable :: precip_regrid(:,:)                      ! Interpolated precip array
  character(len=120) :: fname, name_cmor                    ! Filename variables
  integer           :: ftn
!=== End Variable Definition =======================

  allocate (precip_regrid(LIS_rc%lnc(n),LIS_rc%lnr(n)))
  precip = quad9r
  precip_regrid = quad9r
  fname = name_cmor
  cmorphdata = quad9r
!----------------------------------------------------------------------
! Read raw, uncompressed cmorph data 
!----------------------------------------------------------------------
  ftn = LIS_getNextUnitNumber()
  open(unit=ftn,file=fname, status='old',access='direct', &
      form='unformatted',recl=xd*yd*4,iostat=ios)
  if (ios .eq. 0) then
    read (ftn,rec=1) precip
   
!    
! Set missing values
!
    r=yd
    do i = 1,yd
       do j = 1,xd
         if (precip(j,i) .lt. 0) precip(j,i)=quad9r
       enddo
    enddo
     
!    
! shifting data from 0->360 to -180->+180 LIS standard    
!
    do i = 1,yd
       r=1
       do j = (xd/2)+1,xd
         testout(r,i) = precip(j,i)
         r = r + 1
       enddo
       do j = 1,xd/2
         testout(r,i) = precip(j,i)
         r = r + 1
       enddo
    enddo 
    do i = 1,yd
       do j = 1,xd
         precip(j,i) = testout(j,i)
       enddo
    enddo

!
!=== End of data reconfiguration

  call AGRMET_interp_cmorph(n, xd, yd, precip, LIS_rc%gridDesc(n,:), &
          LIS_rc%lnc(n),LIS_rc%lnr(n),precip_regrid)
     do j = 1,LIS_rc%lnr(n)
        do i = 1,LIS_rc%lnc(n)
           if (precip_regrid(i,j) .ne. -9999.0) then
                 cmorphdata(i,j) = precip_regrid(i,j)   !here is mm/h
           endif
        enddo
     enddo
    
  write(LIS_logunit,*) "Obtained CMORPH precipitation data:: ", fname
 else
   write(LIS_logunit,*) "Missing CMORPH precipitation data ", fname
 endif
 call LIS_releaseUnitNumber(ftn)
 deallocate (precip_regrid)

end subroutine AGRMET_cmorph
