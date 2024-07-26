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
! !ROUTINE: read_nldas2_elev
!  \label{read_nldas2_elev}
!
! !REVISION HISTORY:
!
!  17 Dec 2004; Sujay Kumar; Initial Specificaton
!  24 Aug 2007: Chuck Alonge; Modified for use with NLDAS2 data
!  02 Nov 2012: K. Arsenault; Expanded elev reader for handling
!                          removal of built-in elev corr in NLDAS-2
!
! !INTERFACE:
subroutine read_nldas2_elev( n, findex, narrelev, elevdiff )

! !USES:
  use LDT_coreMod,       only : LDT_rc, LDT_domain
  use LDT_metforcingMod, only : LDT_forc
  use nldas2_forcingMod, only : nldas2_struc
  use LDT_logMod,        only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_fileIOMod,     only : readLISdata

  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
!- Terrain height will be set to run domain:
  real, intent(inout) :: narrelev(LDT_rc%lnc(n),LDT_rc%lnr(n),1)
!- Elevdiff needs to be on forcing domain:
  real, intent(inout) :: elevdiff(nldas2_struc(n)%nc,nldas2_struc(n)%nr,1)

! !DESCRIPTION:
!
!  Opens, reads, and interpolates NLDAS2 model elevation to the LDT
!  grid. The data will be used to perform any topographical 
!  adjustments to the forcing.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[findex]
!   index of the forcing dataset selected
!  \end{description}
! 
!  The routines invoked are: 
!   \begin{description}
!   \end{description}
!EOP
   integer :: ftn1, ftn2
   logical :: file_exists
   integer :: c, r, i
! _____________________________________________________________________________

   narrelev = LDT_rc%udef
   elevdiff = LDT_rc%udef

   if( LDT_rc%lis_map_proj(n).eq."latlon" ) then
    if( LDT_rc%gridDesc(n,4) < (LDT_rc%met_gridDesc(findex,4)-0.0625) .or.  & ! LL Lat
        LDT_rc%gridDesc(n,5) < (LDT_rc%met_gridDesc(findex,5)-0.0625) .or.  & ! LL Lon
        LDT_rc%gridDesc(n,7) > (LDT_rc%met_gridDesc(findex,7)+0.0625) .or.  & ! UR Lat
        LDT_rc%gridDesc(n,8) > (LDT_rc%met_gridDesc(findex,8)+0.0625) ) then  ! UR Lon
      write(LDT_logunit,*)"[ERR] LDT Run domain exceeds NLDAS-2 domain ... ending run."
      call LDT_endrun
    end if
   endif

   if( nldas2_struc(n)%gridDesc(9)  == LDT_rc%gridDesc(n,9) .and. &
       nldas2_struc(n)%gridDesc(10) == LDT_rc%gridDesc(n,10).and. &
       LDT_rc%gridDesc(n,1) == 0 .and. &
       LDT_rc%met_gridtransform_parms(findex) .ne. "neighbor" ) then
      write(LDT_logunit,*) "[WARN]  The NLDAS-2 0.125 deg grid was selected for the"
      write(LDT_logunit,*) "  LDT run domain; however, 'bilinear', 'budget-bilinear',"
      write(LDT_logunit,*) "  or some other unknown option was selected to spatially"
      write(LDT_logunit,*) "  downscale the grid, which will cause errors during runtime."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun()
   endif

   inquire(file = trim(nldas2_struc(n)%file_narrelev), exist=file_exists)
   if(.not. file_exists) then
      write(LDT_logunit,*) "The NLDAS-2/NARR terrain height file ",&
            trim(nldas2_struc(n)%file_narrelev)," is not found."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif
   inquire(file = trim(nldas2_struc(n)%file_elevdiff), exist=file_exists)
   if(.not. file_exists) then
      write(LDT_logunit,*) "NLDAS2 elevation difference file ",&
            trim(nldas2_struc(n)%file_elevdiff)," not found."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif

! -------------------------------------------------------------------
! Open and Read-in Forcing Terrain Hght File - Bring to LIS run domain
! -------------------------------------------------------------------
   write(LDT_logunit,*) "Reading the NLDAS-2/NARR terrain height file: ", &
        trim(nldas2_struc(n)%file_narrelev)

   ftn1 = LDT_getNextUnitNumber()
   open(ftn1, file = nldas2_struc(n)%file_narrelev, form='unformatted',&
        access='direct',recl=4)

   call readLISdata(n, ftn1, LDT_rc%met_proj(findex), &
            LDT_rc%met_gridtransform_parms(findex), LDT_rc%met_gridDesc(findex,:), &
            1, narrelev )

    do r=1,LDT_rc%lnr(n)
       do c=1,LDT_rc%lnc(n)
         if( narrelev(c,r,1) < 0 ) then
            narrelev(c,r,1) = LDT_rc%udef
         endif
      enddo
   enddo
   call LDT_releaseUnitNumber(ftn1)

! -------------------------------------------------------------------
! Open and Read-in Elevation Difference File (to remain on NLDAS2 grid)
! -------------------------------------------------------------------
   write(LDT_logunit,*) "Reading the NLDAS2 elevation difference file: ",&
         trim(nldas2_struc(n)%file_elevdiff)

   ftn2 = LDT_getNextUnitNumber()
   open(ftn2, file = nldas2_struc(n)%file_elevdiff, form='unformatted',&
        access='direct',recl=4)
   i = 0
   do r = 1, nldas2_struc(n)%nr
      do c = 1, nldas2_struc(n)%nc; i = i + 1
         read(ftn2, rec=i) elevdiff(c,r,1)
      enddo
   enddo
   call LDT_releaseUnitNumber(ftn2)


end subroutine read_nldas2_elev
