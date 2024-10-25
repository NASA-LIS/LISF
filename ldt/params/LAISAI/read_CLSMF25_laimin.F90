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
! !ROUTINE: read_CLSMF25_laimin
! \label{read_CLSMF25_laimin}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar;  Initial Specification
!  12 Feb 2013: K. Arsenault; Modified for use with CLSM F2.5 lai parameter
!
! !INTERFACE:
 subroutine read_CLSMF25_laimin(n, array)

! !USES:
  use ESMF
  use LDT_coreMod, only : LDT_rc, LDT_domain
  use LDT_logMod,  only : LDT_logunit, LDT_getNextUnitNumber, &
                          LDT_releaseUnitNumber, LDT_endrun
  use map_utils
  use LDT_laisaiMod

  implicit none

! !ARGUMENTS: 
  integer, intent(in)    :: n    ! nest index
  real,    intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n))

! !DESCRIPTION:
!  This subroutine retrieves the lai fraction climatology for the 
!  specified month and returns the values in the latlon projection
!  
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[array]
!   output field with the retrieved lai fraction
!  \item[maskarray]
!   optional input field of reading in the mask array
!  \end{description}
!
!EOP      
  real, parameter :: min_val = 10000.0
  integer :: ftn
  integer :: c,r,k
  integer :: month
  integer :: glpnr, glpnc, gr, gc
  real    :: param_grid(20)
  real    :: rlat(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real    :: rlon(LDT_rc%lnc(n),LDT_rc%lnr(n))
  logical :: file_exists
  real    :: tmpreal(LDT_rc%nmaskpts(n))
  real, allocatable :: tmparray(:,:)
! __________________________________________________________

  array = min_val
  
  !- Determine global/complete parameter domain number of points:
  param_grid(:) = LDT_rc%mask_gridDesc(n,:)
  glpnr = nint((param_grid(7)-param_grid(4))/param_grid(10)) + 1
  glpnc = nint((param_grid(8)-param_grid(5))/param_grid(9)) + 1
  
  inquire(file=trim(LDT_laisai_struc(n)%laifile), exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) "LAI map ",trim(LDT_laisai_struc(n)%laifile)," not found"
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  endif

  ftn = LDT_getNextUnitNumber()
  open(ftn, file=trim(LDT_laisai_struc(n)%laifile), form="unformatted", status="old" )
  
! Read lai from global file:
  allocate( tmparray(glpnc,glpnr) )
  tmparray = min_val
  
  do month = 1, 12
     read( ftn ) ( tmpreal(k), k=1,LDT_rc%nmaskpts(n) )
     
     k = 0
     do r = 1, glpnr
        do c = 1, glpnc
           if( LDT_rc%global_mask(c,r) > 0. ) then
              k = k + 1
              if( tmpreal(k) .ge. 0 ) then
                 if(tmpreal(k).lt.0.01) tmpreal(k) = 0.01
                 if(tmpreal(k).lt.tmparray(c,r)) then 
                    tmparray(c,r) = tmpreal(k)
                 endif
              endif
           endif
        end do
     enddo
  enddo

  do r=1,glpnr
     do c=1,glpnc
        if(tmparray(c,r).eq.min_val) then 
           tmparray(c,r) = LDT_rc%udef
        endif
     enddo
  enddo

   !- Subset domain:
  do r = 1, LDT_rc%lnr(n)
     do c = 1, LDT_rc%lnc(n)
        call ij_to_latlon(LDT_domain(n)%ldtproj,float(c),float(r),&
             rlat(c,r),rlon(c,r))
        gr = nint((rlat(c,r)-param_grid(4))/param_grid(10))+1
        gc = nint((rlon(c,r)-param_grid(5))/param_grid(9))+1
        array(c,r) = tmparray(gc,gr)
     end do
  end do
  
  deallocate( tmparray)
  
  call LDT_releaseUnitNumber(ftn)

end subroutine read_CLSMF25_laimin
