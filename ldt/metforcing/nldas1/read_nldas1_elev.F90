!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: read_nldas1_elev
!  \label{read_nldas1_elev}
!
! !REVISION HISTORY:
!
!  17 Dec 2004; Sujay Kumar; Initial Specificaton
!  20 Dec 2006; Kristi Arsenault; Changed to read EDAS elevation only
!  03 Dec 2007; Sujay Kumar; Added the abstract method to read the 
!               EDAS data.  
!  08 Aug 2013: K. Arsenault; Expanded elev reader for handling
!                          removal of built-in elev corr in NLDAS-1
!
! !INTERFACE:
subroutine read_nldas1_elev(n, findex, edaselev, elevdiff)

! !USES:
  use LDT_coreMod,       only : LDT_rc, LDT_domain
  use LDT_metforcingMod, only : LDT_forc
  use nldas1_forcingMod, only : nldas1_struc
  use LDT_logMod,        only : LDT_logunit, LDT_getNextUnitNumber, &
          LDT_releaseUnitNumber, LDT_endrun
  use LDT_fileIOMod,     only : readLISdata

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
!- Terrain height will be set to run domain:
  real,  intent(inout) :: edaselev(LDT_rc%lnc(n),LDT_rc%lnr(n),1)
!- Elevdiff needs to be on forcing domain:
  real, intent(inout)  :: elevdiff(nldas1_struc(n)%nc,nldas1_struc(n)%nr,1)

! !DESCRIPTION:
!
!  Opens, reads, and interpolates NLDAS1-EDAS model elevation to the 
!  LDT grid. The data will be used to perform any topographical 
!  adjustments to the forcing.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
! 
!  The routines invoked are: 
!   \begin{description}
!    \item[LDT\_readData](\ref{LDT_readData}) \newline
!     abstract method to read EDAS elevation data in the map 
!     projection used in the LDT model grid. 
!   \end{description}
!EOP

  integer :: c, r
  integer :: i, err
  integer :: nldas1
  integer :: ftn1, ftn2
  logical :: file_exists
! _______________________________________________________________

   edaselev = LDT_rc%udef
   elevdiff = LDT_rc%udef
   nldas1 = nldas1_struc(n)%nc * nldas1_struc(n)%nr

   if( LDT_rc%lis_map_proj.eq."latlon" ) then
    if( LDT_rc%gridDesc(n,4) < LDT_rc%met_gridDesc(findex,4) .or.  & ! LL Lat
        LDT_rc%gridDesc(n,5) < LDT_rc%met_gridDesc(findex,5) .or.  & ! LL Lon
        LDT_rc%gridDesc(n,7) > LDT_rc%met_gridDesc(findex,7) .or.  & ! UR Lat
        LDT_rc%gridDesc(n,8) > LDT_rc%met_gridDesc(findex,8) ) then  ! UR Lon
       write(LDT_logunit,*)" LDT Run domain exceeds NLDAS1 domain ... ending run."
       call LDT_endrun
    end if
   endif

   inquire(file = trim(nldas1_struc(n)%file_edaselev), exist=file_exists)
   if(.not. file_exists) then
      write(LDT_logunit,*) "The NLDAS-1/GTOPO30 terrain height file ",&
            trim(nldas1_struc(n)%file_edaselev)," is not found."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif
   inquire(file = trim(nldas1_struc(n)%file_elevdiff), exist=file_exists)
   if(.not. file_exists) then
      write(LDT_logunit,*) "NLDAS1 elevation difference file ",&
            trim(nldas1_struc(n)%file_elevdiff)," not found."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif

! -------------------------------------------------------------------
!  READ IN Forcing Elevation Difference File (to remain on NLDAS1 grid)
! -------------------------------------------------------------------
   write(LDT_logunit,*) "Reading elevation file ...",&
         trim(nldas1_struc(n)%file_edaselev)

   ftn1 = LDT_getNextUnitNumber()
   open(ftn1,file=nldas1_struc(n)%file_edaselev,form='unformatted', &
        access='direct',recl=4,status='old')

   call readLISdata(n, ftn1, LDT_rc%met_proj(findex), &
            LDT_rc%met_gridtransform_parms(findex), LDT_rc%met_gridDesc(findex,:),&
            1, edaselev)   ! 1 indicates 2D layer

   call LDT_releaseUnitNumber(ftn1)
   write(LDT_logunit,*) "Done reading EDAS elevation file."

! -------------------------------------------------------------------
!- Read in original NLDAS1 EDAS-GTOPO Elevation Difference file::
! -------------------------------------------------------------------
   ftn2 = LDT_getNextUnitNumber()
   open (unit=ftn2, file = nldas1_struc(n)%file_elevdiff, &
         form='unformatted', access='direct', recl=4, iostat=err )

   if( err /= 0 ) THEN
      write(LDT_logunit,*)"STOP: Problem opening elevation difference file"
      write(LDT_logunit,*)"Run without elevation correction option."
      print*, "STOP: problem opening elevation difference file"
      print*, "Try running without elevation correction option."
      call LDT_endrun
   else
      i = 0
      do r = 1, nldas1_struc(n)%nr
         do c = 1, nldas1_struc(n)%nc; i = i + 1
            read(ftn2, rec=i) elevdiff(c,r,1)
         enddo
      enddo
   endif
   call LDT_releaseUnitNumber(ftn2)

   write(LDT_logunit,*) &
        "Read Original NLDAS1 Elevation Difference File:: ", &
         nldas1_struc(n)%file_elevdiff

end subroutine read_nldas1_elev
