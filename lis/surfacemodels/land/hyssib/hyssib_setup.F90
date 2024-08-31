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
! !ROUTINE: hyssib_setup
! \label{hyssib_setup}
!
! !REVISION HISTORY:
!  21 Apr 2004: David Mocko, Conversion from NOAH to HY-SSiB
!  27 Sep 2007: Chuck Alonge, Updates for LIS 5.0
!  27 Oct 2010: David Mocko, changes for HY-SSiB in LIS6.1
!
! !INTERFACE:
subroutine hyssib_setup()
! !USES:
   use LIS_coreMod, only : LIS_rc
   use hyssib_lsmMod
   use hyssibveg_module
   use hyssibalb_module
!
! !DESCRIPTION: 
! 
!  This routine is the entry point to set up the parameters
!  required for Hyssib LSM. These include the vegetation, albedo, 
!  bottom temperature and the initialization of state variables
!  in Hyssib. 
!  
! The routines invoked are: 
! \begin{description}
! \item[hyssibveg\_ini](\ref{hyssibveg_ini}) \newline
!   initializes the Hyssib vegetation parameter module
! \item[hyssibalb\_vegini](\ref{hyssibalb_ini}) \newline
!   initializes the Hyssib albedo parameter module
! \item[hyssib\_setvegparms](\ref{hyssib_setvegparms}) \newline
!   initializes the vegetation-related parameters in Hyssib
! \item[hyssib\_settbot](\ref{hyssib_settbot}) \newline
!   initializes the bottom temperature fields
! \item[hyssib\_settopostd](\ref{hyssib_settopostd}) \newline
!   initializes the std. dev. of topography in hyssib
! \item[hyssib\_gfrac](\ref{hyssib_gfrac}) \newline
!   initializes monthly vegetation parameters in hyssib
! \item[hyssib\_coldstart](\ref{hyssib_coldstart}) \newline
!   initializes the hyssib state variables
! \end{description}
!EOP
   implicit none

   integer :: t, n
!=== End Variable List =================================================

   call hyssibveg_ini()
   call hyssibalb_ini()
   call hyssib_setvegparms()
   call hyssib_settbot()
   call hyssib_settopostd()
   call hyssib_gfrac()

   call hyssib_coldstart()
   
   do n=1,LIS_rc%nnest
      hyssib_struc(n)%count = 0 
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
         hyssib_struc(n)%hyssib(t)%green = 0
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
         hyssib_struc(n)%hyssib(t)%uwind = 0
         hyssib_struc(n)%hyssib(t)%vwind = 0
         hyssib_struc(n)%hyssib(t)%rainf_in = 0
         hyssib_struc(n)%hyssib(t)%snowf_in = 0
         hyssib_struc(n)%hyssib(t)%tair = 0
         hyssib_struc(n)%hyssib(t)%qair = 0
         hyssib_struc(n)%hyssib(t)%psurf = 0
         hyssib_struc(n)%hyssib(t)%swdown = 0
         hyssib_struc(n)%hyssib(t)%lwdown = 0         
      enddo
   enddo

end subroutine hyssib_setup

