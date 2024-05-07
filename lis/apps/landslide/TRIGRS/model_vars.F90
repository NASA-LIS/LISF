!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

! By Rex L. Baum, 1 April 2004
module model_vars
	integer,parameter:: double=kind(1d0)
	integer:: nts,kper,nmax1,nmax2,nmn,nmin
	integer,allocatable:: jsav(:), ix(:),jy(:) ! added ix() & jy() 4/21/2010 
	real:: test1,dg2rad
	real,allocatable:: q(:),qb(:)
	real (double):: eps,tmin,tmax,ts,tinc,qt,tns,beta,qmax !
	real (double):: test,nodat,sumex,dusz,dcf,vf0,p0zmx ! Added p0zmx 6 May 2013, RLB
	real (double):: nodat_Slope ! SY: For assisting rain assignment from 2-d LIS grid to 1-d TRIGRS array
	real (double):: celsiz,param(6),parami(6),ti,tis,pi,smt,lard,xllc,yllc,zmn(1),zmx(1) ! added xllc & yllc 12/24/2010
	real (double),allocatable:: p(:),ptran(:),pzero(:),bline(:),chi(:) ! added chi 12/22/2010 RLB
	real (double),allocatable:: r(:),fc(:),fw(:),thz(:),kz(:),tcap(:)
	real (double),allocatable:: trz(:),uwsp(:),gs(:),qtime(:),qts(:) ! Added qts 21 Feb 2013, RLB 
end module model_vars

