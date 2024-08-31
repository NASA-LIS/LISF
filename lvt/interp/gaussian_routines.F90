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
! !ROUTINE: gaussian_comp_lats
! \label{gaussian_comp_lats}
!
! !INTERFACE:
subroutine gaussian_comp_lats(jmax, gaussian_lat_array)
! !USES:

   implicit none
! !ARGUMENTS:
   integer, intent(in)                :: jmax
   real, dimension(jmax), intent(out) :: gaussian_lat_array

!
! !DESCRIPTION: 
!  This routine computes the latitudes of the global quasi-regular Gaussian 
!  grid corresponding to the total number of latitude circles given by jmax.
!
!   \begin{description}
!   \item [jmax]
!      the total number of latitude circles
!   \item [gaussian\_lat\_array]
!      array of computed latitudes, in degrees
!   \end{description}
!EOP

   real*8, parameter :: pi=3.14159265358979

   real, allocatable, dimension(:) :: slat, wlat

   real*8  :: dpr
   integer :: j

   dpr =180./pi

   allocate(slat(jmax))
   allocate(wlat(jmax))

   call gausslat(jmax,slat,wlat)

   do j = 1, jmax
      gaussian_lat_array(j) = dpr*asin(slat(jmax-j+1))
   enddo

   deallocate(slat)
   deallocate(wlat)

end subroutine gaussian_comp_lats


!BOP
! !ROUTINE: gaussian_find_row
! \label{gaussian_find_row}
!
! !INTERFACE:
function gaussian_find_row(jmax, gaussian_lat_array, lat)
! !USES:

   implicit none
! !ARGUMENTS:
   integer, intent(in)               :: jmax
   real, dimension(jmax), intent(in) :: gaussian_lat_array
   real, intent(in)                  :: lat
!
! !DESCRIPTION: 
!  This function computes the row number within the global quasi-regular 
!  Gaussian grid corresponding to the latitude given by lat.
!
!   \begin{description}
!   \item [jmax]
!      the total number of latitude circles
!   \item [gaussian\_lat\_array]
!      array of computed latitudes, in degrees
!   \item [lat]
!      latitude to search for
!   \end{description}
!EOP

   integer :: gaussian_find_row

   real    :: eps
   integer :: r

   eps = 180./(2*jmax)

   gaussian_find_row = -9999
   do r = 1, jmax
      if ( abs(gaussian_lat_array(r) - lat) < eps ) then
      !if ( gaussian_latitudes(1,r) == lat ) then
         gaussian_find_row = r
      endif
   enddo

   if ( gaussian_find_row == -9999 ) then
      write(*,*) '[ERR] gaussian_find_row -- '// &
                               'Could not find lat',lat
      stop
   endif

end function gaussian_find_row
