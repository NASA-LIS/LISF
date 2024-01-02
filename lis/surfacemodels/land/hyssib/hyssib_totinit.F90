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
! !ROUTINE: hyssib_totinit
!  \label{hyssib_totinit}
!
! !REVISION HISTORY:
!  14 Jun 2002: Sujay Kumar, Initial Specification
!  21 Apr 2004: David Mocko, Conversion from NOAH to HY-SSiB
!  28 Aug 2007: Chuck Alonge, Updates for LIS 5.0
!  27 Oct 2010: David Mocko, changes for HY-SSiB in LIS6.1
!
! !INTERFACE:
subroutine hyssib_totinit(n)
! !USES:
   use hyssib_lsmMod   
   use LIS_coreMod, only : LIS_rc

   implicit none

   integer, intent(in) :: n
!
! !DESCRIPTION:
!
!  This routine resets the time-averaged Hyssib variables after a 
!  model output. The routine is called after a model output call. 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP

   integer t
!=== End Variable List =================================================

   do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
      hyssib_struc(n)%hyssib(t)%swnet = 0
      hyssib_struc(n)%hyssib(t)%lwnet = 0
      hyssib_struc(n)%hyssib(t)%qle = 0
      hyssib_struc(n)%hyssib(t)%qh = 0
      hyssib_struc(n)%hyssib(t)%qg = 0
      hyssib_struc(n)%hyssib(t)%qf = 0
      hyssib_struc(n)%hyssib(t)%qv = 0
      hyssib_struc(n)%hyssib(t)%qtau = 0
      hyssib_struc(n)%hyssib(t)%qa = 0
      hyssib_struc(n)%hyssib(t)%delsurfheat = 0
      hyssib_struc(n)%hyssib(t)%delcoldcont = 0
      hyssib_struc(n)%hyssib(t)%snowf = 0
      hyssib_struc(n)%hyssib(t)%rainf = 0
      hyssib_struc(n)%hyssib(t)%evap = 0
      hyssib_struc(n)%hyssib(t)%qs = 0
      hyssib_struc(n)%hyssib(t)%qrec = 0
      hyssib_struc(n)%hyssib(t)%qsb = 0
      hyssib_struc(n)%hyssib(t)%qsm = 0
      hyssib_struc(n)%hyssib(t)%qfz = 0
      hyssib_struc(n)%hyssib(t)%qst = 0
      hyssib_struc(n)%hyssib(t)%delsoilmoist = 0
      hyssib_struc(n)%hyssib(t)%delswe = 0
      hyssib_struc(n)%hyssib(t)%delintercept = 0
      hyssib_struc(n)%hyssib(t)%snowt = 0
      hyssib_struc(n)%hyssib(t)%vegtc = 0
      hyssib_struc(n)%hyssib(t)%baresoilt = 0
      hyssib_struc(n)%hyssib(t)%avgsurft = 0
      hyssib_struc(n)%hyssib(t)%radteff = 0
      hyssib_struc(n)%hyssib(t)%albedo = 0
      hyssib_struc(n)%hyssib(t)%swe = 0
      hyssib_struc(n)%hyssib(t)%sweveg = 0
      hyssib_struc(n)%hyssib(t)%soilmoist(1) = 0
      hyssib_struc(n)%hyssib(t)%soilmoist(2) = 0
      hyssib_struc(n)%hyssib(t)%soilmoist(3) = 0
      hyssib_struc(n)%hyssib(t)%soilmoist1m = 0
      hyssib_struc(n)%hyssib(t)%soiltemp = 0
      hyssib_struc(n)%hyssib(t)%soilwet = 0
      hyssib_struc(n)%hyssib(t)%soilwetrz = 0
      hyssib_struc(n)%hyssib(t)%potevap = 0
      hyssib_struc(n)%hyssib(t)%ecanop = 0
      hyssib_struc(n)%hyssib(t)%tveg = 0
      hyssib_struc(n)%hyssib(t)%esoil = 0
      hyssib_struc(n)%hyssib(t)%rootmoist = 0
      hyssib_struc(n)%hyssib(t)%canopint = 0
      hyssib_struc(n)%hyssib(t)%subsnow = 0
      hyssib_struc(n)%hyssib(t)%subsurf = 0
      hyssib_struc(n)%hyssib(t)%acond = 0
      hyssib_struc(n)%hyssib(t)%ccond = 0
      hyssib_struc(n)%hyssib(t)%snowfrac = 0
      hyssib_struc(n)%hyssib(t)%snowdepth = 0
      hyssib_struc(n)%hyssib(t)%sliqfrac = 0
      hyssib_struc(n)%hyssib(t)%ect = 0
      hyssib_struc(n)%hyssib(t)%eci = 0
      hyssib_struc(n)%hyssib(t)%egs = 0
      hyssib_struc(n)%hyssib(t)%hc = 0
      hyssib_struc(n)%hyssib(t)%hg = 0
      hyssib_struc(n)%hyssib(t)%radnvisdir = 0
      hyssib_struc(n)%hyssib(t)%radnvisdif = 0
      hyssib_struc(n)%hyssib(t)%radnnirdir = 0
      hyssib_struc(n)%hyssib(t)%radnnirdif = 0
      hyssib_struc(n)%hyssib(t)%salbvisdir = 0
      hyssib_struc(n)%hyssib(t)%salbvisdif = 0
      hyssib_struc(n)%hyssib(t)%salbnirdir = 0
      hyssib_struc(n)%hyssib(t)%salbnirdif = 0
      hyssib_struc(n)%hyssib(t)%radtc = 0
      hyssib_struc(n)%hyssib(t)%radtg = 0
      hyssib_struc(n)%hyssib(t)%watcan = 0
      hyssib_struc(n)%hyssib(t)%watgrd = 0
      hyssib_struc(n)%hyssib(t)%snocan = 0
      hyssib_struc(n)%hyssib(t)%snogrd = 0
      hyssib_struc(n)%hyssib(t)%wet1 = 0
      hyssib_struc(n)%hyssib(t)%wet2 = 0
      hyssib_struc(n)%hyssib(t)%wet3 = 0
      hyssib_struc(n)%hyssib(t)%swdown_out = 0
      hyssib_struc(n)%hyssib(t)%lwdown_out = 0
      hyssib_struc(n)%hyssib(t)%snowtcount = 0
      hyssib_struc(n)%hyssib(t)%albedocount = 0
      hyssib_struc(n)%hyssib(t)%sliqfraccount = 0
   enddo
   hyssib_struc(n)%count = 0
end subroutine hyssib_totinit

