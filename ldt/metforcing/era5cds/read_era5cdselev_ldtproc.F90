!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: read_era5cdselev_ldtproc
! \label{read_era5cdselev_ldtproc}
!
! !REVISION HISTORY:
!
!  17Dec2004; Sujay Kumar; Initial Specificaton
!
! !INTERFACE:
subroutine read_era5cdselev_ldtproc(n, findex, change)
! !USES:
  use LDT_constantsMod,   only : LDT_CONST_PATH_LEN
  use LDT_coreMod,        only : LDT_rc, LDT_domain
  use LDT_metforcingMod,  only : LDT_forc
  use LDT_fileIOMod,      only : LDT_read_param
  use LDT_logMod,         only : LDT_logunit, LDT_endrun
  use era5cds_forcingMod, only : era5cds_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
  integer, intent(in) :: change

! !DESCRIPTION:
!
!  Opens, reads, and interpolates ERA5CDS model elevation to the LDT
!  grid. The data will be used to perform any topographical 
!  adjustments to the forcing.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP
  integer :: c,r,line1,line2,nc_dom,line
  integer :: glnc, glnr
  real    :: go(LDT_rc%lnc(n),LDT_rc%lnr(n))

  if ( trim(LDT_rc%met_ecor(findex)).eq."lapse-rate")then 
     write(LDT_logunit,*) " Reading the ERA5CDS elevation: "
     if ( change == 0 ) then ! constant for all period
        call LDT_read_param(n,"ELEV_ERA5CDS",go)
     else
        write(LDT_logunit,*) 'ERR: invalid update request. ', change, &
                         'is not for ERA5CDS elev'
        call LDT_endrun
     endif

     do r=1,LDT_rc%lnr(n)
        do c=1,LDT_rc%lnc(n)
           if(LDT_domain(n)%gindex(c,r).ne.-1) then 
              LDT_forc(n,findex)%modelelev(LDT_domain(n)%gindex(c,r)) = go(c,r)
           endif
        enddo
     enddo

  endif

end subroutine read_era5cdselev_ldtproc
