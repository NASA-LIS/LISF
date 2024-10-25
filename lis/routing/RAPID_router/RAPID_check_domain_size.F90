!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !DESCRIPTION: This subroutine checks the size of static data for RAPID.
!              
!  Reference:
!
! !REVISION HISTORY:
! Mar 10, 2021: Yeosang Yoon, Initial Implementation
!
! !USES:
subroutine RAPID_check_domain_size(n)

  use RAPID_routingMod

  implicit none

  integer, intent(in)     :: n

  integer                 :: eof, status
  integer                 :: m 
  character(len=200)      :: buf

  ! check total number of river connect information
  open(10,file=trim(RAPID_routing_struc(n)%connectfile),status='old',action='read')
  m=1;
  do
     read(10,*,iostat=eof) buf
     if (eof/=0) exit

     m=m+1
  end do
  m=m-1 ! adjust size
  close(10)
  RAPID_routing_struc(n)%n_riv_tot=m
 
  ! check total number of river basins
  open(11,file=trim(RAPID_routing_struc(n)%basinIDfile),status='old',action='read')
  m=1;
  do
     read(11,*,iostat=eof) buf
     if (eof/=0) exit

     m=m+1
  end do
  m=m-1 ! adjust size
  close(11)
  RAPID_routing_struc(n)%n_riv_bas=m

  ! check number of reach in weight table file
  open(45,file=trim(RAPID_routing_struc(n)%weightfile),status='old',action='read')
  read(45,'(A)',iostat=eof) buf  ! read header in weight table
  m=1;
  do
     read(45,*,iostat=eof) buf
     if (eof/=0) exit

     m=m+1
  end do
  m=m-1 ! adjust size
  close(45)
  RAPID_routing_struc(n)%n_wei_table=m

end subroutine RAPID_check_domain_size
