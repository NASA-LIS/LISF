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
!
! !ROUTINE: read_TRMM3B42V6.F90
! \label{read_TRMM3B42V6}
!
! !REVISION HISTORY: 
!  17 Jul 2001: Jon Gottschalck; Initial code
!  04 Feb 2002: Jon Gottschalck; Added necessary code to use global precip
!               observations with LDT_domain 3 (2x2.5)
!  30 Jul 2002: Jon Gottschalck; Added code to use 3B42V6 and Persiann precip
!               data
!  25 Aug 2006: Yudong Tian; Modified for LDT 4.2 release
!  19 Apr 2013: Jonathan Case; Corrected some initializations and comparison
!               against undefined value
!  21 Jun 2013: Soni Yatheendradas; changes from earlier code to avoid
!               (a) alternate file skip,
!               (b) jump to previous day TRMM, and
!               (c) absence of rain rate weighting
!
! !INTERFACE:
subroutine read_TRMM3B42V6 (n, fname, findex, order, ferror_TRMM3B42V6)
! !USES:
  use LDT_coreMod, only           : LDT_rc, LDT_domain
  use LDT_logMod, only            : LDT_logunit, LDT_getNextUnitNumber, &
                                    LDT_releaseUnitNumber
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_metforcingMod, only     : LDT_forc
  use TRMM3B42V6_forcingMod, only : TRMM3B42V6_struc

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  character(len=LDT_CONST_PATH_LEN) :: fname
  integer, intent(in) :: findex
  integer, intent(in) :: order
  integer             :: ferror_TRMM3B42V6
  !integer             :: filehr ! SY

!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  TRMM 3B42V6 data and interpolates to the LDT domain.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[fname]
!    name of the 3-hourly TRMM 3B42V6 file
!  \item[ferror\_TRMM3B42V6]
!    flag to indicate success of the call (=0 indicates success)
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!  \item[interp\_TRMM3B42V6](\ref{interp_TRMM3B42V6}) \newline
!    spatially interpolates the TRMM 3B42V6 data
!  \end{description}
!EOP

  integer :: index1, nTRMM3B42V6

  !==== Local Variables=======================

  integer :: ios
  integer :: i,j,xd,yd
  parameter(xd=1440,yd=400)                            ! Dimension of original 3B42V6 data

  real :: precip(xd,yd)                                ! Original real precipitation array
  logical*1,allocatable  :: lb(:)
  real, allocatable :: precip_regrid(:,:)                  ! Interpolated precipitation array
  integer :: ftn

  !=== End Variable Definition =======================

  !------------------------------------------------------------------------
  ! Fill necessary arrays to assure not using old 3B42V6 data
  !------------------------------------------------------------------------
  ! J.Case (4/22/2013) -- Make consistent with Stg4/NMQ routines
  if(order.eq.1) then
     LDT_forc(n,findex)%metdata1 = LDT_rc%udef ! J.Case
  elseif(order.eq.2) then 
     LDT_forc(n,findex)%metdata2 = LDT_rc%udef ! J.Case
  endif
  allocate (precip_regrid(LDT_rc%lnc(n),LDT_rc%lnr(n)))
!  precip = -1.0
!  if(order.eq.1) then
!     LDT_forc(n,findex)%metdata1 = -1.0
!  elseif(order.eq.2) then 
!     LDT_forc(n,findex)%metdata2 = -1.0
!  endif
  precip_regrid = -1.0 ! J.Case
  !------------------------------------------------------------------------
  ! Find 3B42V6 precip data, read it in and assign to forcing precip array.
  ! Must reverse grid in latitude dimension to be consistent with LDAS grid
  !------------------------------------------------------------------------
  ftn = LDT_getNextUnitNumber()
  open(unit=ftn,file=fname, status='old', &
       &          access='direct',recl=xd*yd*4, &
       &          form='unformatted',iostat=ios)

  if (ios .eq. 0) then
     read (ftn,rec=1) precip
     Do j=1, xd
        Do i=1, yd
           if (precip(j, i) .LT. 0.0 ) precip(j, i) = 0.0    ! reset to 0 for weird values
        End Do
     End Do

! J.Case (4/19/2013) -- Test print out of raw precip array
! write (99,*) precip

     !------------------------------------------------------------------------
     ! Interpolating to desired LDT_domain and resolution
     ! Global precip datasets not used currently to force NLDAS
     !------------------------------------------------------------------------
     !print*, "Writing un-interpolated 3B42V6 precipitation out "
     !open(71, file="TRMM3B42V6-ungrid.1gd4r", access="direct", &
     !    recl=xd*yd*4, form="unformatted")
     ! write(71, rec=1) precip
     !close(71)

     nTRMM3B42V6 = TRMM3B42V6_struc(n)%nc*TRMM3B42V6_struc(n)%nr
     allocate(lb(nTRMM3B42V6))
     lb = .true.
     call interp_TRMM3B42V6(n, nTRMM3B42V6, precip, lb, LDT_rc%gridDesc, &
          !LDT_rc%lnc(n),LDT_rc%lnr(n),precip_regrid) ! SY
          LDT_rc%lnc(n),LDT_rc%lnr(n),precip_regrid, findex) ! SY
     deallocate (lb) 

     !print*, "Writing interpolated 3B42V6 precipitation out "
     !open(73, file="TRMM3B42V6-regrid.1gd4r", access="direct", &
     !    recl=LDT_rc%d%lnr*LDT_rc%d%lnc*4, form="unformatted")
     ! write(73, rec=1) precip_regrid
     !close(73)
     !print*, "Writing interpolated 3B42V6 precipitation out finished"

! J.Case (4/19/2013) -- Test print out of the regridded precip (on LDT grid).
! write (98,*) precip_regrid

     do j = 1,LDT_rc%lnr(n)
        do i = 1,LDT_rc%lnc(n)
           if (precip_regrid(i,j) .ne. -1.0) then
              index1 = LDT_domain(n)%gindex(i,j)
              if(index1 .ne. -1) then
                 if(order.eq.1) then 
                    LDT_forc(n,findex)%metdata1(1,index1) = precip_regrid(i,j)   !here is mm/h
                 elseif(order.eq.2) then 
                    LDT_forc(n,findex)%metdata2(1,index1) = precip_regrid(i,j)   !here is mm/h
                 endif
              endif
           endif
        enddo
     enddo

! J.Case (4/19/2013) -- Test print out of the suppdata precip (on LDT grid).
! write (97,*) LDT_forc(n,findex)%metdata1(1,:)

     ferror_TRMM3B42V6 = 1
     write(LDT_logunit,*) "Obtained 3B42 V6 precipitation data ", trim(fname)
  else
     write(LDT_logunit,*) "Missing 3B42 V6 precipitation data ", trim(fname)
     ferror_TRMM3B42V6 = 0
  endif
  call LDT_releaseUnitNumber(ftn)

  deallocate (precip_regrid)

end subroutine read_TRMM3B42V6
