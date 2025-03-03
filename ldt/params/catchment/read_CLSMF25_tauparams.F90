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
! !ROUTINE: read_CLSMF25_tauparams
!  \label{read_CLSMF25_tauparams}

! !REVISION HISTORY:
!  25 Nov 2012: K. Arsenault; Initial Specification
!
! !INTERFACE:
subroutine read_CLSMF25_tauparams(n,array1,array2,maskarray) 

! !USES:
  use LDT_coreMod,       only : LDT_rc, LDT_domain
  use LDT_logMod,        only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use map_utils
  use CLSMF25_parmsMod
!EOP      
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  real, intent(inout) :: array1(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real, intent(inout) :: array2(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real, intent(inout) :: maskarray(LDT_rc%lnc(n),LDT_rc%lnr(n))
!
! !DESCRIPTION:
!  This subroutine retrieves the Catchment LSM surface layer timescale
!  (tau) parameters from a text-based formatted file and places it on
!  a separate spatially-distributed data field.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[array1]
!   output field for atau parameter
!  \item[array2]
!   output field for btau parameter
!  \item[maskarray]
!   input mask field for mapping parameter to Catchment map
!  \end{description}
!
!EOP      
  integer  :: ftn
  logical  :: file_exists
  integer  :: read_status
  integer  :: c, r, k, m
  real     :: dzsf
  integer  :: glpnr, glpnc, gr, gc
  real     :: param_grid(20)
  real     :: rlat(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real     :: rlon(LDT_rc%lnc(n),LDT_rc%lnr(n))

  integer  :: tmptileid, dummy_int
  real     :: tmpreal(4)
  real, allocatable :: tmparray1(:,:)
  real, allocatable :: tmparray2(:,:)
! _____________________________________________________________________

   array1 = LDT_rc%udef
   array2 = LDT_rc%udef

 ! NOTE: Units of dz** are [mm] while excess/deficits from catchment()
 !       are in SI units (ie kg/m^2) or loosely speaking, in mm of water.
 !       In other words, density of water (1000 kg/m^3) is built 
 !       into dz** (reichle, 5 Feb 04).
 !- Surface layer:
   dzsf = CLSMF25_struc(n)%dzsfcrd * 1000. 
 ! ----

!- Determine global/complete parameter domain number of points:
   param_grid(:) = LDT_rc%mask_gridDesc(n,:)
   glpnr = nint((param_grid(7)-param_grid(4))/param_grid(10)) + 1
   glpnc = nint((param_grid(8)-param_grid(5))/param_grid(9)) + 1

   inquire(file=trim(CLSMF25_struc(n)%sltsfile), exist=file_exists)
   if(.not.file_exists) then 
     write(LDT_logunit,*) "Catchment F2.5 surface layer timescale  &
               (tau) parameter file ",trim(CLSMF25_struc(n)%sltsfile)," not found."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
   endif
  
   ftn = LDT_getNextUnitNumber()
   write(LDT_logunit,*) "[INFO] Reading surface layer timescale (tau) data &
             from: ",trim(CLSMF25_struc(n)%sltsfile)
   open(ftn, file=CLSMF25_struc(n)%sltsfile, form="formatted")

!- Loop over the CLSM mask values and read in each line from the tau file:
   allocate( tmparray1(glpnc,glpnr) )
   allocate( tmparray2(glpnc,glpnr) )
   tmparray1 = LDT_rc%udef
   tmparray2 = LDT_rc%udef

! - For future subsetted domains:
!   do r = 1, LDT_rc%lnr(n)
!      do c = 1, LDT_rc%lnc(n)
!         if( maskarray(c,r) > 0. ) then
!- For now - complete domains:
   do r = 1, glpnr
      do c = 1, glpnc
         if( LDT_rc%global_mask(c,r) > 0. ) then

         !- Surface layer timescale parameters:
            read(ftn,*,iostat=read_status) tmptileid, &
                 dummy_int, (tmpreal(m), m=1,4)

         ! Select atau and btau depending on surface layer depth:

           if( abs(dzsf-20.) < 1e-4 ) then   ! use atau2, btau2
!              array1(c,r) = tmpreal(1)   
!              array2(c,r) = tmpreal(2)   
              tmparray1(c,r) = tmpreal(1)   
              tmparray2(c,r) = tmpreal(2)   

           elseif( abs(dzsf-50.) < 1e-4 ) then ! use atau5, btau5
!              array1(c,r) = tmpreal(3)   
!              array2(c,r) = tmpreal(4)   
              tmparray1(c,r) = tmpreal(3)   
              tmparray2(c,r) = tmpreal(4)   
           else
             write (LDT_logunit,*) "[ERR] read_CLSMF25_tauparms : "
             write (LDT_logunit,*) "   unknown value for dzsf " 
             call LDT_endrun()
           end if

         end if
      end do
   end do 

!- Subset domain:
   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n)
         call ij_to_latlon(LDT_domain(n)%ldtproj,float(c),float(r),&
                           rlat(c,r),rlon(c,r))
         gr = nint((rlat(c,r)-param_grid(4))/param_grid(10))+1
         gc = nint((rlon(c,r)-param_grid(5))/param_grid(9))+1
         array1(c,r) = tmparray1(gc,gr)
         array2(c,r) = tmparray2(gc,gr)
      end do
   end do
   deallocate( tmparray1 )
   deallocate( tmparray2 )


  call LDT_releaseUnitNumber(ftn)

end subroutine read_CLSMF25_tauparams
