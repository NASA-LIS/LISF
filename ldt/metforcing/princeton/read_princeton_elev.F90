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
! !ROUTINE: read_princeton_elev
! \label{read_princeton_elev}
!
! !REVISION HISTORY:
!
!  1 Feb 2007; Hiroko Kato; Initial Specificaton 
! 25 May 2014; KR Arsenault/H. Kato-Beaudoing, added Princeton elev field
!
! !INTERFACE:
subroutine read_princeton_elev(n, findex, princetonelev, elevdiff)

! !USES:
  use LDT_coreMod,         only : LDT_rc, LDT_domain
  use LDT_metforcingMod,   only : LDT_forc
  use LDT_logMod,          only : LDT_logunit, LDT_getNextUnitNumber, &
           LDT_releaseUnitNumber, LDT_endrun
  use princeton_forcingMod,only : princeton_struc
  use LDT_fileIOMod,       only : LDT_transform_paramgrid

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
!- Terrain height will be set to run domain:
  real, intent(inout) :: princetonelev(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real, intent(inout) :: elevdiff(LDT_rc%met_nc(findex), LDT_rc%met_nr(findex))

!
! !DESCRIPTION:
!
!  Opens and reads PRINCETON model elevation to the LDT
!  grid. The data will be used to perform any topographical 
!  adjustments to the forcing.
!  The elevation file needs to be preprocessed to fit the running 
!  resolution and domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!
!  The routines invoked are: 
!   \begin{description}
!   \item[LDT\_readData](\ref{LDT_readData}) \newline
!     call the abstract method to get the elevation data in the
!     LDT projection. 
!   \end{description}
!EOP

  integer   :: ftn
  integer   :: i, j, jj
  integer   :: c, r
  integer   :: inpts, outpts
  logical   :: file_exists
  integer,parameter :: nx=360, ny=180
  real      :: read_data(nx,ny), yrev_data(nx,ny)

  real      :: elev1d(LDT_rc%met_nc(findex)*LDT_rc%met_nr(findex))
  logical*1 :: lb(LDT_rc%met_nc(findex)*LDT_rc%met_nr(findex))

  real      :: elev_regrid(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  logical*1 :: lb_regrid(LDT_rc%lnc(n)*LDT_rc%lnr(n))

! ___________________________________________________

  princetonelev = LDT_rc%udef
  elevdiff = LDT_rc%udef

!  if ( trim(LDT_rc%met_ecor(findex)) .ne."none") then 

     inquire(file=princeton_struc(n)%elevfile, exist=file_exists)
     if(.not. file_exists) then
        write(LDT_logunit,*) "The PRINCETON terrain height file, ",&
              trim(princeton_struc(n)%elevfile),", is not found."
        write(LDT_logunit,*) "Program stopping ..."
        call LDT_endrun
     endif
     write(LDT_logunit,*) "Reading the PRINCETON elevation map: ", &
                           trim(princeton_struc(n)%elevfile)
     
   ! Read in free format ascii data:
     ftn = LDT_getNextUnitNumber()
     open(ftn,file=princeton_struc(n)%elevfile,form='formatted',status='old')
     do j=1, 180
        read(ftn,*) read_data(:,j)
     enddo
     call LDT_releaseUnitNumber(ftn)
 
!     print*,minval(read_data),maxval(read_data)
!     write(*,*) '120:',read_data(:,120)
!     write(*,*) '180:',read_data(:,180)
 
   ! First grid at north-west corner (180W,90N)
   !  Need to flip north-south direction to be in LIS format
     do j=1, 180
        jj = 180 - j + 1
        yrev_data(:,j) = read_data(:,jj)
     enddo

  !- Initialize arrays for grid transformation:
     inpts  = LDT_rc%met_nc(findex)*LDT_rc%met_nr(findex)
     outpts = LDT_rc%lnr(n)*LDT_rc%lnc(n)
     elev1d = 0.0
     lb     = .false.
     lb_regrid = .true.

  !- Assign 2-D array to 1-D for aggregation routines:
     i = 0
     do r = 1, LDT_rc%met_nr(findex)
        do c = 1, LDT_rc%met_nc(findex);  i = i + 1
           elev1d(i) = yrev_data(c,r)
           if( elev1d(i) .ne. LDT_rc%udef ) lb(i) = .true.
        enddo
     enddo

  !- Interp elevation field to output field:
     call LDT_transform_paramgrid(n, LDT_rc%met_gridtransform_parms(findex), &
              LDT_rc%met_gridDesc(findex,:), inpts, 1, elev1d, lb, &
              outpts, elev_regrid, lb_regrid )

  !- Convert 1D to 2D elevation file:
     i = 0
     do r = 1, LDT_rc%lnr(n)
        do c = 1, LDT_rc%lnc(n)
           i = i + 1
           princetonelev(c,r) = elev_regrid(i)
        end do
     end do

!  endif

end subroutine read_princeton_elev
