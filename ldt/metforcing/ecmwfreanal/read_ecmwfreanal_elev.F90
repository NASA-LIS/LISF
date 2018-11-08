!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: read_ecmwfreanal_elev
!  \label{read_ecmwfreanal_elev}
!
! !REVISION HISTORY:
!  17Dec2004; Sujay Kumar; Initial Specificaton
!  3 Dec2007; Sujay Kumar; Added the abstract method to read the 
!               EDAS data.  
!
! !INTERFACE:
subroutine read_ecmwfreanal_elev(n, findex, elev, elevdiff)

! !USES:
  use LDT_coreMod,           only : LDT_rc, LDT_domain
  use LDT_metforcingMod,     only : LDT_forc
  use LDT_logMod,            only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_fileIOMod,         only : readLISdata
  use ecmwfreanal_forcingMod, only : ecmwfreanal_struc

!EOP
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
!- Terrain height will be set to run domain:
  real, intent(inout) :: elev(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real, intent(inout) :: elevdiff(LDT_rc%met_nc(findex),LDT_rc%met_nr(findex))
!
! !DESCRIPTION:
!  Opens, reads, and interpolates ECMWF reanalysis model elevation. 
!  The data will be used to perform any topographical adjustments 
!  to the forcing. 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!
!  The routines invoked are: 
!   \begin{description}
!    \item[readLISdata](\ref{readLISdata}) \newline
!     abstract method to read ECMWM-Reanal elevation data  
!     in the map projection used in the LDT model grid. 
!   \end{description}
!EOP
  integer  :: c,r
  integer  :: ftn
  logical  :: file_exists
  real     :: gridDesci(20)
  character(50) :: proj
  real, allocatable :: go(:,:,:)
! _______________________________________________

  elev = LDT_rc%udef
  elevdiff = LDT_rc%udef

  inquire(file=ecmwfreanal_struc(n)%elevfile, exist=file_exists)
  if(.not. file_exists) then
     write(LDT_logunit,*) "The ECMWF-Reanalysis terrain height file, ",&
           trim(ecmwfreanal_struc(n)%elevfile),", is not found."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  endif
  write(LDT_logunit,*) 'Reading the ECMWF Reanalysis elevation ',&
       ecmwfreanal_struc(n)%elevfile

  select case( ecmwfreanal_struc(n)%elevtransform )

   case( "neighbor", "average", "none" )
     write(LDT_logunit,*) " ECMWF Reanalysis Elevation File Grid, being transformed with:"
     write(LDT_logunit,*) trim(ecmwfreanal_struc(n)%elevtransform)
   case default
     write(LDT_logunit,*) " ERR: Since the version of ECMWF Reanalysis Elevation "
     write(LDT_logunit,*) "  file being read-in has a resolution at 0.01 deg (globally),"
     write(LDT_logunit,*) "  the user should select either 'neighbor' or 'average' if going"
     write(LDT_logunit,*) "  to transform or aggregate the grid.  'none' also works if staying"
     write(LDT_logunit,*) "  at 0.01 deg lat-lon grid."
     write(LDT_logunit,*) " Current option selected: ",trim(ecmwfreanal_struc(n)%elevtransform)
     write(LDT_logunit,*) " Program stopping ..."
     call LDT_endrun
   end select

   gridDesci = 0
   gridDesci(1) = 0. 
   gridDesci(2) = 36000
   gridDesci(3) = 15000
   gridDesci(4) = -59.995000
   gridDesci(5) = -179.995000
   gridDesci(6) = 128
   gridDesci(7) = 89.99500
   gridDesci(8) = 179.99500
   gridDesci(9) = 0.0100
   gridDesci(10) =  0.0100
   gridDesci(20) = 64

   ftn = LDT_getNextUnitNumber()
   open(ftn,file=ecmwfreanal_struc(n)%elevfile,form='unformatted',&
        access='direct',recl=4)

   allocate( go(LDT_rc%lnc(n),LDT_rc%lnr(n),1) )
   go = LDT_rc%udef
   proj = "latlon"
   call readLISdata(n, ftn, proj, ecmwfreanal_struc(n)%elevtransform,&
                    gridDesci(:), 1, go )
  
   call LDT_releaseUnitNumber(ftn)

   elev(:,:) = go(:,:,1)
   deallocate(go)


end subroutine read_ecmwfreanal_elev
