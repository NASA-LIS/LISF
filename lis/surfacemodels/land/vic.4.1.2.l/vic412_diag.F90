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
! !ROUTINE: vic412_diagnoseoutputvar
! \label{vic412_diagnoseoutputvar}
!
! !REVISION HISTORY:
! Nov 2011; James Geiger, Initial implementation of VIC 4.1.1 into LIS.
! 
! !INTERFACE:
subroutine vic412_diagnoseoutputvar(nest, tile, moc_index, vlevel, value, &
                                    units, lenu, direction, lend)
! !USES:
   use LIS_coreMod, only : LIS_rc
   use LIS_histDataMod

   implicit none
! !ARGUMENTS: 
   integer, intent(in)             :: nest, tile, moc_index, vlevel
   integer, intent(in)             :: lenu, lend
   real*8, intent(in)              :: value
   character(len=lenu), intent(in) :: units
   character(len=lend), intent(in) :: direction
!
! !DESCRIPTION:
! This routine calls LIS\_diagnoseOutputVar to process VIC's output variables.
!
!  The arguments are: 
!  \begin{description}
!  \item[nest]
!   index of the nest
!  \item[tile]
!   index of the tile
!  \item[moc\_index]
!   index identifying which output variable to process
!  \item[vlevel]
!   vertical level to process
!  \item[value]
!   output value to process
!  \item[units]
!   units of the output variable
!  \item[lenu]
!   length of the units string
!  \item[direction]
!   direction of the output variable
!  \item[lend]
!   length of the direction string
!  \end{description}
!EOP

   ! call LIS_diagnoseOutputVar(nest, tile, moc_index, value=real(value), &
   !                           vlevel=vlevel, unit=units, direction=direction)
   call LIS_diagnoseSurfaceOutputVar(nest,tile,moc_index,value=real(value),   &
                              vlevel=vlevel, unit=units, direction=direction, &
                              surface_type=LIS_rc%lsm_index)

end subroutine vic412_diagnoseoutputvar
