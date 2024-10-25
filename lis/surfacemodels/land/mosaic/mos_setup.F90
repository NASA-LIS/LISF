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
! !ROUTINE: mos_setup
! \label{mos_setup}
!
! !REVISION HISTORY:
! 3 Jun  2003: Jon Gottschalck; Initial Code
! 25 Sep 2007: Sujay Kumar; Upgraded for LIS5.0
!
! !INTERFACE:
subroutine mos_setup()
! !USES:
  use LIS_coreMod, only : LIS_rc
  use mos_lsmMod

  implicit none
!
! !DESCRIPTION: 
! 
!  This routine is the entry point to set up the parameters
!  required for Mosaic LSM. These include the soils and topography 
!  and the initialization of state variables
!  in Mosaic. 
!  
! The routines invoked are: 
! \begin{description}
! \item[setmosp](\ref{setmosp}) \newline
!   initializes the static soil and topography fields 
! \item[mapsib2umd](\ref{mapsib2umd}) \newline
!   maps the SIB vegetation types to UMD classification
! \item[mos\_coldstart](\ref{mos_coldstart}) \newline
!   initializes the noah state variables
! \end{description}
!EOP

  integer :: t, n

  do n=1,LIS_rc%nnest
     call setmosp(n)
     call mapsib2umd()
     call mos_coldstart(n)
   
     do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        mos_struc(n)%mos(t)%swnet = 0
        mos_struc(n)%mos(t)%lwnet = 0
        mos_struc(n)%mos(t)%qle = 0
        mos_struc(n)%mos(t)%qh = 0
        mos_struc(n)%mos(t)%qg = 0
        mos_struc(n)%mos(t)%rainf = 0
        mos_struc(n)%mos(t)%snowf = 0
        mos_struc(n)%mos(t)%evap = 0
        mos_struc(n)%mos(t)%qs = 0
        mos_struc(n)%mos(t)%qsm = 0
        mos_struc(n)%mos(t)%qsb = 0
        mos_struc(n)%mos(t)%swe = 0
        mos_struc(n)%mos(t)%snowcover = 0 
        mos_struc(n)%mos(t)%soilmoist1 = 0
        mos_struc(n)%mos(t)%soilmoist2 = 0
        mos_struc(n)%mos(t)%soilmoist3 = 0
        mos_struc(n)%mos(t)%soilwet = 0
        mos_struc(n)%mos(t)%ecanop = 0
        mos_struc(n)%mos(t)%tveg = 0
        mos_struc(n)%mos(t)%esoil = 0
        mos_struc(n)%mos(t)%rootmoist = 0
        mos_struc(n)%mos(t)%canopint = 0
        mos_struc(n)%mos(t)%soilm_prev = 0
        mos_struc(n)%mos(t)%swe_prev = 0
        mos_struc(n)%mos(t)%dtcanal = 0 
        mos_struc(n)%count = 0
     enddo
  enddo
  
 end subroutine mos_setup

