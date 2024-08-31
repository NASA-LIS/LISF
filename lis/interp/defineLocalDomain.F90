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
! !ROUTINE: defineLocalDomain
!  \label{defineLocalDomain}
! 
! !REVISION HISTORY: 
! 
! 10 July 2007: Sujay Kumar; Initial Specification
! 
! !INTERFACE: 
subroutine defineLocalDomain(lisGD, inGD, outGD)
! !USES: 
  use LIS_logMod,  only : LIS_logunit, LIS_endrun
  
  implicit none

  real,      intent(IN)    :: lisGD(50)
  real,      intent(IN)    :: inGD(50)
  real,      intent(OUT)   :: outGD(50)
!
! !DESCRIPTION: 
!  This routine generates a description of the sub-domain based on the 
!  specifications of the running domain, native resolution and the 
!  output domain. Derives the lat/lon extents from the lis grid description,
!  the map projection and resolution to be used from the input/native 
!  grid description. 
!
!EOP
  real                    :: lat1,lon1,lat2,lon2
  real                    :: lat1_in, lon1_in, lat2_in, lon2_in
  integer                 :: lnc, lnr
  
  outGD = 0 

  if(lisGD(1).eq.0) then 

     outGD(1) = inGD(1)

     if(inGD(1).eq.0) then ! from lat/lon to a lat/lon domain
        lat1 = lisGD(4)
        lon1 = lisGD(5)

        lat2 = lisGD(7)
        lon2 = lisGD(8)
        
! find these lat/lons in the input GD. 
        
        lat1_in = floor((lat1-inGD(3))/inGD(9))*inGD(9)+inGD(3)
        lon1_in = floor((lon1-inGD(4))/inGD(10))*inGD(10)+inGD(4)
        
!        print*, lat1, 'becomes ',lat1_in
!        print*, lon1, 'becomes ',lon1_in

        lat2_in = ceiling((lat2-inGD(3))/inGD(9))*inGD(9)+inGD(3)
        lon2_in = ceiling((lon2-inGD(4))/inGD(10))*inGD(10)+inGD(4)

!        print*, lat2, 'becomes ',lat2_in
!        print*, lon2, 'becomes ',lon2_in
       
        lnc = nint((lat2_in-lat1_in)/inGD(9))+1
        lnr = nint((lon2_in-lon1_in)/inGD(10))+1

        outGD(2)  = lnc
        outGD(3)  = lnr
        outGD(4)  = lat1_in
        outGD(5)  = lon1_in
        outGD(6)  = 128
        outGD(7)  = lat2_in
        outGD(8)  = lon2_in
        outGD(9)  = inGD(9)
        outGD(10) = inGD(10)
!        print*, 'here ',outGD(1:20)
!        stop
     else
        write(LIS_logunit,*) 'the use of subsetting in this projection is '
        write(LIS_logunit,*) 'not currently supported ... Program stopping....'
     endif
  else
     write(LIS_logunit,*) 'the use of subsetting in this projection is '
     write(LIS_logunit,*) 'not currently supported ... Program stopping....'
     call LIS_endrun()
  endif

end subroutine defineLocalDomain
