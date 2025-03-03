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
! !ROUTINE: read_CLSMF25_lai
! \label{read_CLSMF25_lai}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!  20 Feb 2006: Sujay Kumar; Modified to support nesting
!  12 Feb 2013: KR Arsenault; Modified for Catchment LSM LAI
!  19 Jul 2021: Eric Kemp; Removed third argument (not used)
!
! !INTERFACE:
!subroutine read_CLSMF25_lai(n, array, maskarray)
subroutine read_CLSMF25_lai(n, array)

! !USES:
  use ESMF
  use LDT_coreMod, only : LDT_rc, LDT_domain
  use LDT_logMod,  only : LDT_logunit, LDT_getNextUnitNumber, &
                          LDT_releaseUnitNumber, LDT_endrun
  use LDT_laisaiMod
  use map_utils

  implicit none

! !ARGUMENTS:
  integer, intent(in)    :: n
  real,    intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n),12)
  !real, optional, intent(inout) :: maskarray(LDT_rc%lnc(n),LDT_rc%lnr(n))

! !DESCRIPTION:
!  This subroutine retrieves the leaf area index (LAI) climatology for the
!  specified month and returns the values for Catchment F2.5 LSM.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   time index (month or quarter)
!  \item[array]
!   output field with the retrieved LAI
!  \end{description}
!
!EOP
  integer :: ftn
  integer :: c, r, k
  integer :: month
  integer :: glpnr, glpnc, gr, gc
  real    :: param_grid(20)
  real    :: rlat(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real    :: rlon(LDT_rc%lnc(n),LDT_rc%lnr(n))
  logical :: file_exists
  real    :: tmpreal(LDT_rc%nmaskpts(n))
  real, allocatable :: tmparray(:,:,:)
! ________________________________________________________

   array = LDT_rc%udef

!- Determine global/complete parameter domain number of points:
   param_grid(:) = LDT_rc%mask_gridDesc(n,:)
   glpnr = nint((param_grid(7)-param_grid(4))/param_grid(10)) + 1
   glpnc = nint((param_grid(8)-param_grid(5))/param_grid(9)) + 1

   inquire(file=trim(LDT_laisai_struc(n)%laifile), exist=file_exists)
   if (.not. file_exists) then
      write(LDT_logunit,*) "[ERR] LAI map ", &
           trim(LDT_laisai_struc(n)%laifile)," not found"
      write(LDT_logunit,*) "[ERR] Program stopping ..."
      call LDT_endrun
   endif

   ftn = LDT_getNextUnitNumber()
   open(ftn, file=trim(LDT_laisai_struc(n)%laifile), &
        form="unformatted",status='old' )

   !- Read LAI from global file:
   allocate( tmparray(glpnc,glpnr,12) )
   tmparray = LDT_rc%udef

   do month = 1, 12
      read( ftn ) ( tmpreal(k), k=1, LDT_rc%nmaskpts(n) )

      k = 0
! - For future subsetted domains:
!      do r = 1, LDT_rc%lnr(n)
!         do c = 1, LDT_rc%lnc(n)
!            if( maskarray(c,r) > 0. ) then
! - For now - complete domains:
      do r = 1, glpnr
         do c = 1, glpnc
            if ( LDT_rc%global_mask(c,r) > 0. ) then

              k = k + 1
              if( tmpreal(k) < 0 ) then
                 tmparray(c,r,month) = LDT_rc%udef
              else
                 tmparray(c,r,month) = tmpreal(k)
              endif
           endif
           if ( tmparray(c,r,month) < 0.01 ) tmparray(c,r,month) = 0.01

         end do
      enddo

   !- Subset domain:
      do r = 1, LDT_rc%lnr(n)
         do c = 1, LDT_rc%lnc(n)
            call ij_to_latlon(LDT_domain(n)%ldtproj, float(c), float(r),&
                              rlat(c,r), rlon(c,r))
            gr = nint((rlat(c,r)-param_grid(4))/param_grid(10))+1
            gc = nint((rlon(c,r)-param_grid(5))/param_grid(9))+1
            array(c,r,month) = tmparray(gc,gr,month)
         end do
      end do

   end do  ! end month loop
   deallocate( tmparray )

   call LDT_releaseUnitNumber(ftn)

end subroutine read_CLSMF25_lai
