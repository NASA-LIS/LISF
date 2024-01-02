!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "vic411_atmos_forcing.h"
!BOP
! !ROUTINE: vic411_read_gridded_forcing_data
! \label{vic411_read_gridded_forcing_data}
!
! !INTERFACE:
subroutine vic411_read_gridded_forcing_data(n, findex, filename, ferror)
! !USES:
   use LIS_coreMod,        only : LIS_rc, LIS_domain
   use LIS_logMod,         only : LIS_getNextUnitNumber, &
                                  LIS_releaseUnitNumber, &
                                  LIS_logunit, LIS_endrun
   use LIS_metforcingMod,  only : LIS_forc
   use vic_forcingMod,     only : vicforcing_struc
   use vic411_lsmMod,      only : vic411_struc

! !ARGUMENTS: 
   implicit none

   integer, intent(in)            :: n
   integer, intent(in)            :: findex
   character(len=*), intent(in) :: filename
   integer, intent(out)           :: ferror

! !DESCRIPTION: 
! This routine reads the VIC-processed forcing data from disk and
! stores it within the suppdata arrays.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[findex]
!   index of the supplemental forcing scheme
!  \end{description}
!EOP

   integer :: NC, NR
   real, allocatable, dimension(:) :: forcing_array, tempvic

   integer :: ftn, i, j, c, r, index

   logical :: exists
   
   inquire(file=trim(filename), exist=exists)
   if ( exists ) then
      NC  = vicforcing_struc(n)%NC
      NR  = vicforcing_struc(n)%NR

      allocate(forcing_array(NC*NR))
      allocate(tempvic(LIS_rc%ngrid(n)))

      write(LIS_logunit, *) "VIC forcing: opening ", trim(filename)

      ftn = LIS_getNextUnitNumber()
      open(ftn,file=trim(filename),form='unformatted',access='direct',&
           recl=4*NC*NR)

      do j = 1, NUM_ATMOS_FORCING
         read(ftn, rec=j) forcing_array

         call interp_vic(n, findex, j, NC*NR, forcing_array, tempvic)

         if ( j == ATMOS_AIR_TEMP+1) then
            if ( vic411_struc(n)%debugging_convert_units == 1 ) then
            ! convert to ALMA units
               where ( tempvic /= -9999.0 )
                  tempvic = tempvic + VIC_KELVIN
               endwhere
            endif
         endif

         if ( j == ATMOS_PREC+1) then
            if ( vic411_struc(n)%debugging_convert_units == 1 ) then
            ! convert to ALMA units
               where ( tempvic /= -9999.0 )
                  tempvic = tempvic / (vicforcing_struc(n)%forcingInterval*3600)
               endwhere
            endif
         endif

         do r = 1, LIS_rc%lnr(n)
            do c = 1, LIS_rc%lnc(n)
               index = LIS_domain(n)%gindex(c,r)
               if ( index .ne. -1 ) then 
                  vicforcing_struc(n)%metdata1(j,index) = tempvic(index)
               endif
            enddo
         enddo

      enddo

      deallocate(forcing_array)
      deallocate(tempvic)

      close(ftn)
      call LIS_releaseUnitNumber(ftn)

      ferror = 0
   else
      ferror = 1
   endif

end subroutine vic411_read_gridded_forcing_data

!BOP
! 
! !ROUTINE: interp_vic
! \label{interp_vic}
! 
! !INTERFACE: 
subroutine interp_vic(n,findex,var,NLEN,f,tempvic)
! !USES: 
   use LIS_coreMod,     only : LIS_rc, LIS_domain
   use vic_forcingMod,  only : vicforcing_struc

   implicit none
! !ARGUMENTS: 
   integer, intent(in)                           :: n 
   integer, intent(in)                           :: findex
   integer, intent(in)                           :: var
   integer, intent(in)                           :: NLEN 
   real, dimension(NLEN), intent(in)             :: f
   real, dimension(LIS_rc%ngrid(n)), intent(out) :: tempvic
! 
! !DESCRIPTION: 
!  Interpolates the VIC formatted forcing data onto the LIS grid. 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[NLEN]
!    size of the VIC forcing data
!  \item[var]
!    variable index
!  \item[f]
!    native VIC forcing data 
!  \item[g]
!    interpolated VIC forcing data
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
!  \end{description}
!EOP
   real                      :: gridDesco(50)    ! Input,output grid info arrays
   integer                   :: count1,mo,iret,i,j
   logical*1                 :: lb(NLEN)
   logical*1                 :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
   real                      :: go(LIS_rc%lnc(n)*LIS_rc%lnr(n))

!------------------------------------------------------------------------     
! Initializing input and output grid arrays
!------------------------------------------------------------------------
   gridDesco = 0
   gridDesco = LIS_rc%gridDesc(n,:)
   mo = LIS_rc%lnc(n) * LIS_rc%lnr(n)
!------------------------------------------------------------------------
! Defining input data bitmap
!------------------------------------------------------------------------
   lb = .true.
   do i = 1, NLEN
     if ( f(i) == LIS_rc%udef ) then
         lb(i) = .false.
     endif
   enddo
!------------------------------------------------------------------------
! Defining output data bitmap
!------------------------------------------------------------------------
   lo = .true.
!------------------------------------------------------------------------
! Interpolate data from VIC forcing grid to LIS grid
!------------------------------------------------------------------------
   if(LIS_rc%met_interp(findex).eq."bilinear") then 
      call bilinear_interp(gridDesco,lb,f,lo,go,NLEN,mo,         &
           LIS_domain(n)%lat, LIS_domain(n)%lon,&
           vicforcing_struc(n)%w111,vicforcing_struc(n)%w121,            &
           vicforcing_struc(n)%w211,vicforcing_struc(n)%w221,            &
           vicforcing_struc(n)%n111,vicforcing_struc(n)%n121,            &
           vicforcing_struc(n)%n211,vicforcing_struc(n)%n221,LIS_rc%udef,&
           iret)
   elseif(LIS_rc%met_interp(findex).eq."budget-bilinear") then 
      if ( var == ATMOS_PREC+1) then
         call conserv_interp(gridDesco,lb,f,lo,go,NLEN,mo,          &
              LIS_domain(n)%lat, LIS_domain(n)%lon,&
              vicforcing_struc(n)%w112,vicforcing_struc(n)%w122,            &
              vicforcing_struc(n)%w212,vicforcing_struc(n)%w222,            &
              vicforcing_struc(n)%n112,vicforcing_struc(n)%n122,            &
              vicforcing_struc(n)%n212,vicforcing_struc(n)%n222,LIS_rc%udef,&
              iret)
      else
         call bilinear_interp(gridDesco,lb,f,lo,go,NLEN,mo,         &
              LIS_domain(n)%lat, LIS_domain(n)%lon,&
              vicforcing_struc(n)%w111,vicforcing_struc(n)%w121,            &
              vicforcing_struc(n)%w211,vicforcing_struc(n)%w221,            &
              vicforcing_struc(n)%n111,vicforcing_struc(n)%n121,            &
              vicforcing_struc(n)%n211,vicforcing_struc(n)%n221,LIS_rc%udef,&
              iret)
      endif
   endif

   count1 = 0
   do j = 1, LIS_rc%lnr(n)
      do i = 1, LIS_rc%lnc(n)
         if( LIS_domain(n)%gindex(i,j) /= -1 ) then
            tempvic(LIS_domain(n)%gindex(i,j)) = go(i+count1)
         endif
      enddo
      count1 = count1 + LIS_rc%lnc(n)
   enddo

end subroutine interp_vic
