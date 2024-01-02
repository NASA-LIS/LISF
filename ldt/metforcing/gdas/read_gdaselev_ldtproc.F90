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
! !ROUTINE: read_gdaselev_ldtproc
! \label{read_gdaselev_ldtproc}
!
! !REVISION HISTORY:
!
!  17Dec2004; Sujay Kumar; Initial Specificaton
!  17Feb2015; KR Arsenault; Adapted for LDT
!
! !INTERFACE:
subroutine read_gdaselev_ldtproc(n, findex, change)
! !USES:
  use LDT_coreMod,       only : LDT_rc, LDT_domain
  use LDT_metforcingMod, only : LDT_forc
  use LDT_fileIOMod,     only : LDT_read_param
  use LDT_logMod,        only : LDT_logunit, LDT_endrun
  use gdas_forcingMod,   only : gdas_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
  integer, intent(in) :: change
! !DESCRIPTION:
!
!  Opens, reads, and interpolates GDAS model elevation to the LDT
!  grid. The data will be used to perform any topographical 
!  adjustments to the forcing.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[LDT\_readData](\ref{LDT_readData}) \newline
!    Abstract method to read the elevation of the forcing
!    data in the same map projection used in LDT.
!  \end{description}
!EOP
  integer :: c,r
  real :: go(LDT_rc%lnc(n),LDT_rc%lnr(n))

  if ( LDT_rc%met_ecor(findex).eq."lapse-rate") then 
     write(LDT_logunit,*) 'Reading the GDAS elevation'
     if ( change == 0 ) then ! period 1980--1991
        ! Note that for this time period we have a difference file
        ! instead of a proper elevation file.  Correct this.
        call LDT_read_param(n,"ELEV_GDAS_T126",go)
     elseif ( change == 1 ) then ! period 1991--2000
        call LDT_read_param(n,"ELEV_GDAS_T126",go)
     elseif ( change == 2 ) then ! period 2000--2002
        call LDT_read_param(n,"ELEV_GDAS_T170",go)
     elseif ( change == 3 ) then ! period 2002--2005
        call LDT_read_param(n,"ELEV_GDAS_T254",go)
     elseif ( change == 4 ) then ! period 2005--2010
        call LDT_read_param(n,"ELEV_GDAS_T382",go)
     elseif ( change == 5 ) then ! period 2010--2015
        call LDT_read_param(n,"ELEV_GDAS_T574",go)
     elseif ( change == 6 ) then ! period 2015--
        call LDT_read_param(n,"ELEV_GDAS_T1534",go)
     else
        write(LDT_logunit,*) 'ERR: invalid update request. ', change, &
                         'is not in {1,2,3,4,5,6}'
        call LDT_endrun
     endif

     if ( change == 0 ) then ! period 1980--1991
        ! Note that for this time period we have a difference file
        ! instead of a proper elevation file.  Correct this.
        !
        ! For now set LDT_domain%grid%elev = 0, and set the modelelev to the
        ! additive inverse of the difference data.
        ! We want to compute x - y, but we are reading (x-y) instead of y.
        go = -go
     endif

     do r=1,LDT_rc%lnr(n)
        do c=1,LDT_rc%lnc(n)
           if(LDT_domain(n)%gindex(c,r).ne.-1) then 
              LDT_forc(n,findex)%modelelev(LDT_domain(n)%gindex(c,r)) = go(c,r)
           endif
        enddo
     enddo

  endif

end subroutine read_gdaselev_ldtproc
