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
! !ROUTINE: read_CLSMF25_albvis
!  \label{read_CLSMF25_albvis}

! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar;  Initial Specification
!  12 Feb 2013: KR Arsenault; Modified for CLSM albedo scale factors
!
! !INTERFACE:
 subroutine read_CLSMF25_albvis(n,albvr,albvf,maskarray) 

! !USES:
  use LDT_coreMod, only : LDT_rc, LDT_domain
  use LDT_logMod,  only : LDT_logunit, LDT_getNextUnitNumber, &
                          LDT_releaseUnitNumber, LDT_endrun
  use map_utils
  use CLSMF25_parmsMod

  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n
  real, intent(inout) :: albvr(LDT_rc%lnc(n),LDT_rc%lnr(n),12)
  real, intent(inout) :: albvf(LDT_rc%lnc(n),LDT_rc%lnr(n),12)
  real, optional, intent(inout) :: maskarray(LDT_rc%lnc(n),LDT_rc%lnr(n))
!
! !DESCRIPTION:
!  This subroutine retrieves the albedo VIS factors for CLSM F2.5 model.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[albvr]
!   output field with the retrieved albedo VIS direct factors
!  \item[albvf]
!   output field with the retrieved albedo VIS diffuse factors
!  \end{description}
!
!EOP      
  integer  :: ftn1
  logical  :: file_exists
  integer  :: c, r, t
  integer  :: month
  integer  :: glpnr, glpnc, gr, gc
  real     :: param_grid(20)
  real     :: rlat(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real     :: rlon(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real     :: tmpreal(LDT_rc%nmaskpts(n))

  real, allocatable :: tmparray1(:,:,:)
  real, allocatable :: tmparray2(:,:,:)
! ___________________________________________________________

   albvr = LDT_rc%udef
   albvf = LDT_rc%udef

!- Determine global/complete parameter domain number of points:
   param_grid(:) = LDT_rc%mask_gridDesc(n,:)
   glpnr = nint((param_grid(7)-param_grid(4))/param_grid(10)) + 1
   glpnc = nint((param_grid(8)-param_grid(5))/param_grid(9)) + 1

   inquire(file=trim(CLSMF25_struc(n)%albvisfile), exist=file_exists)
   if(.not.file_exists) then
      write(LDT_logunit,*) "Albedo VIS file ",trim(CLSMF25_struc(n)%albvisfile)," not found."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif
!   select case ( LDT_rc%alb_gridtransform(n) )
!     case( "none" )
!       write(LDT_logunit,*) "[INFO] Reading CLSM F2.5 Albedo VIS "
!     case default
!       write(LDT_logunit,*) "[ERR] The spatial transform option selected for CLSM F2.5 Albedo VIS"
!       write(LDT_logunit,*) "     monthly is not recognized nor recommended.  Please select:  none"
!       write(LDT_logunit,*) "Program stopping ..."
!       call LDT_endrun
!   end select
  
   ftn1 = LDT_getNextUnitNumber()
   open(ftn1, file=CLSMF25_struc(n)%albvisfile,form='unformatted',status='old' )

!- Read CLSM albedo VIS scale factor values:  
   allocate( tmparray1(glpnc,glpnr,12), tmparray2(glpnc,glpnr,12) )
   tmparray1 = LDT_rc%udef
   tmparray2 = LDT_rc%udef

   do month = 1, 12

      read (ftn1) (tmpreal(t), t=1,LDT_rc%nmaskpts(n))

      t = 0
! - For future subsetted domains:
!      do r = 1, LDT_rc%lnr(n)
!        do c = 1, LDT_rc%lnc(n)
!           if( maskarray(c,r) > 0 ) then
! - For now - complete domains:
      do r = 1, glpnr
         do c = 1, glpnc
            if( LDT_rc%global_mask(c,r) > 0. ) then

              t = t + 1
              if( tmpreal(t) < 0 ) then
!                 albvr(c,r,month) = LDT_rc%udef
!                 albvf(c,r,month) = LDT_rc%udef
                 tmparray1(c,r,month) = LDT_rc%udef
                 tmparray2(c,r,month) = LDT_rc%udef
              else
!                 albvr(c,r,month) = tmpreal(t)
!                 albvf(c,r,month) = tmpreal(t)
                 tmparray1(c,r,month) = tmpreal(t)
                 tmparray2(c,r,month) = tmpreal(t)
              endif
           endif
!           if( albvr(c,r,month) < 0.01 ) albvr(c,r,month) = 0.01
!           if( albvf(c,r,month) < 0.01 ) albvf(c,r,month) = 0.01
           if( tmparray1(c,r,month) < 0.01 ) tmparray1(c,r,month) = 0.01
           if( tmparray2(c,r,month) < 0.01 ) tmparray2(c,r,month) = 0.01

        end do
      enddo

   !- Subset domain:
      do r = 1, LDT_rc%lnr(n)
         do c = 1, LDT_rc%lnc(n)
            call ij_to_latlon(LDT_domain(n)%ldtproj,float(c),float(r),&
                             rlat(c,r),rlon(c,r))
            gr = nint((rlat(c,r)-param_grid(4))/param_grid(10))+1
            gc = nint((rlon(c,r)-param_grid(5))/param_grid(9))+1
            albvr(c,r,month) = tmparray1(gc,gr,month)
            albvf(c,r,month) = tmparray2(gc,gr,month)
         end do
      end do

   end do  ! End month loop
   deallocate( tmparray1, tmparray2 )

   call LDT_releaseUnitNumber(ftn1)

end subroutine read_CLSMF25_albvis
