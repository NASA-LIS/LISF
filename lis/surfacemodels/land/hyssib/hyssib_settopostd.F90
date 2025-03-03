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
! !ROUTINE: hyssib_settopostd
! \label{hyssib_settopostd}
!
! !REVISION HISTORY:
! 03 Jun 2005: David Mocko, Conversion from NOAH to HY-SSiB
! 27 Sep 2007: Chuck ALonge, Updates for LIS 5.0
!
! !INTERFACE:
subroutine hyssib_settopostd
! !USES:
  use LIS_coreMod, only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_logunit
  use hyssib_lsmMod           
!
! !DESCRIPTION:
!  This subroutine retrieves standard deviation of topography data
!  for the Hyssib LSM.  Currently static data is being used. 
!EOP
  implicit none

!=== Local Variables ===================================================

  integer :: n,i,c,r !loop counters
  integer :: line,line1,line2,glnc,glnr
  real, allocatable :: placetopostd(:,:)

  do n=1,LIS_rc%nnest

!-----------------------------------------------------------------------
! Read in std. dev. topography fields
!-----------------------------------------------------------------------
  allocate(placetopostd(LIS_rc%lnc(n),LIS_rc%lnr(n)))

  write(LIS_logunit,*) 'Opening Hyssib TOPOSTD File: ',trim(hyssib_struc(n)%topostdfile)

  open(12,file=hyssib_struc(n)%topostdfile,status='old',&
       access='direct',form="unformatted",recl=4)
  !! currently will work only for lat/lon projection - need to add as a 
  !! parameter plugin in LIS in order to handle other projections 
  line1 = nint((LIS_rc%gridDesc(n,4)-LIS_rc%gridDesc(n,34))/LIS_rc%gridDesc(n,10))+1
  line2 = nint((LIS_rc%gridDesc(n,5)-LIS_rc%gridDesc(n,35))/LIS_rc%gridDesc(n,9))+1
  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        glnc = line2+c-1
        glnr = line1+r-1
        line = (glnr-1)*nint(LIS_rc%gridDesc(n,32))+glnc
        read(12,rec=line) placetopostd(c,r)
     enddo
  enddo
  close(12)
  
  write(LIS_logunit, * ) 'MSG: sethyssibp -- Read TOPOSTD file'
  do i = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if (placetopostd(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row).ne.LIS_rc%udef)&
          then
        hyssib_struc(n)%hyssib(i)%tempstd = &
             placetopostd(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col,&
             LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     endif
  enddo
  
  deallocate(placetopostd)

  enddo

end subroutine hyssib_settopostd

