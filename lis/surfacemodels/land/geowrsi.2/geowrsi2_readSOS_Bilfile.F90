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
! !ROUTINE:  geowrsi2_readSOS_Bilfile
!  \label{geowrsi2_readSOS_Bilfile}
!
! !REVISION HISTORY:
! 25 Apr 2014: KR Arsenault;  Implemented separate SOS-file reader
!
! !INTERFACE:
 subroutine geowrsi2_readSOS_Bilfile( n, numpts, year, SOS, SOSa, SOS_offset )

! !USES:
  use ESMF
  use LIS_coreMod,  only : LIS_rc, LIS_domain, &
                           LIS_masterproc, LIS_npes
  use LIS_logMod,   only : LIS_logunit, LIS_endrun
  use geowrsi2_lsmMod, only: geowrsi2_struc
  use geowrsi2_physics_module,   only : offsetTStepByNumTSteps
  use geowrsi2_arraymgmt_module, only : nullify_ptr, alloc_arr, dealloc_arr
  use fbil_module

  implicit none

! !ARGUMENTS: 
  integer, intent(in)      :: n     ! nest
  integer, intent(in)      :: numpts
  integer, intent(in)      :: year
  integer*4, intent(inout) :: SOS(numpts)
  integer*4, intent(inout) :: SOSa(numpts)
  integer*4, intent(in)    :: SOS_offset
  
!
! !DESCRIPTION:
! 
!  This routine reads the SOS and SOS anomaly fields from BIL
!  formatted files.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!EOP
  integer             :: t, c, r, l, rc
  character(len=4)    :: nest_str
  character(len=4)    :: cyr4
  type(charN)         :: sos_filename
  type(charN)         :: sosa_filename
  character(300),target :: temp_filename

! BIL parameters/variables:
  logical             :: fileread_ok
  logical             :: real4ptr2dIsNull
  real*4, pointer, dimension(:,:) :: real4ptr2d

! ________________________________________________________

!- Setup BIL Geographic Coordination information (from header files):
   allocate(gCoords)

 ! Set the grid space
   call geowrsi2_set_gcoords(n)
   call calcXyBounds( &
                      addOffset_arg=(0 +1),     &
                      ifMidPointRoundUp=.true., &
                      minLon=gCoords%minLon,    &
                      ulxmap=gCoords%minLon,    &
                      xdim=gCoords%pixLon,      &
                      maxLon=gCoords%maxLon,    &
                      ulymap=gCoords%maxLat,    &
                      minLat=gCoords%minLat,    &
                      ydim=gCoords%pixLat,      &
                      maxLat=gCoords%maxLat,    &
                      minX=gCoords%minX,        &
                      maxX=gCoords%maxX,        &
                      minY=gCoords%minY,        &
                      maxY=gCoords%maxY         )

 ! Allocate temporary memory for staging the data in:
   call nullify_ptr(real4ptr2dIsNull, real4ptr2d=real4ptr2d)
   real4ptr2dIsNull = alloc_arr(real4ptr2dIsNull, &
                                dim1Sz=((gCoords%maxX-gCoords%minX) + 1), &
                                dim2Sz=((gCoords%maxY-gCoords%minY) + 1), &
                                real4ptr2d=real4ptr2d)

 ! Initialize SOS fields:
   SOS  = LIS_rc%udef     ! 0 dekad init
   SOSa = 0               ! 0 set by B. Wind

 ! Account for nest domain in SOS/SOSa filenames:
   write(unit=nest_str, fmt='(a2,i2.2)') '.d',n
   write(unit=cyr4,fmt='(i4.4)')  year

!== SOS Timestep (in dekad) Distributed field  

 ! Read SOS-calculated file(s) (based on growing season):
   temp_filename = trim(LIS_rc%odir)//"/SOS_current_"//cyr4//nest_str
   sos_filename%str => temp_filename
   fileread_ok = populate_arr_from_FEWSNET_bil_2d( &
                         mapFilename=sos_filename, &   
                         real4ptr2d=real4ptr2d )
   if( fileread_ok .eqv. .false. )  real4ptr2d = LIS_rc%udef
   if( fileread_ok ) then
      write(LIS_logunit,*) "WRSI-Parameter: SOS current file (",&
                trim(temp_filename),") read-in successfully."
!      print *, "WRSI-Parameter: SOS current file (",&
!                trim(temp_filename),") read-in successfully."
   else
      write(LIS_logunit,*) "ERR MSG:  SOS current NOT read-in (file: ",&
                           trim(temp_filename),")."
      write(LIS_logunit,*) "End run being called ..."
      call LIS_endrun
   endif

 ! "Add or subtract a user-specified number of [time steps] to/from the SOS
 !  at each pixel.  This is useful if [one wishes] to increment SOS over a period
 !  of time and see the variations in water balance given that planting was done at
 !  SOS, SOS+1, SOS+2 etc."  From the GeoWRSI in-program 'Help' dialog.

   do t = 1, numpts
      SOS(t) = int(real4ptr2d(LIS_domain(n)%tile(t)%col,LIS_domain(n)%tile(t)%row))

   !- Apply SOS_offset to SOS field in WRSI run-mode:
      if((SOS_offset == int(LIS_rc%udef)) .or. &
         (SOS(t) == int(LIS_rc%udef))          &
        ) cycle

      if( (SOS_offset /= LIS_rc%udef) .and. &
          (SOS_offset /= 0)           .and. &
          (SOS(t) /= LIS_rc%udef)     .and. &
          (SOS(t) >= 1)               .and. &
          (SOS(t) <= 36)              .and. &
        ! Another check: Don't let the SOS become less than the initial
        ! timestep - don't subtract from InitialTStep
          (.not. (SOS(t) == geowrsi2_struc(n)%wrsi(t)%InitialTStep .and. SOS_offset < 0) ) &
         ) then

        ! Tweak SOS by offset amount to protect against SOS being >36 or <1
          call offsetTStepByNumTSteps( SOS(t), SOS_offset )
      endif
   end do


!== SOS Anomaly Distributed Field (based on SOS and SOSClim):

   temp_filename = trim(LIS_rc%odir)//"/SOS_anomaly_"//cyr4//nest_str
   sosa_filename%str => temp_filename
   fileread_ok = populate_arr_from_FEWSNET_bil_2d( &
                    mapFilename=sosa_filename, &   ! from lis.config file
                    real4ptr2d=real4ptr2d )
   if( fileread_ok .eqv. .false. )  real4ptr2d = LIS_rc%udef
   do t = 1, LIS_rc%npatch(n,LIS_rc%lsm_index)
      SOSa(t) = int(real4ptr2d(LIS_domain(n)%tile(t)%col,LIS_domain(n)%tile(t)%row))
   end do
   if( fileread_ok ) then
      write(LIS_logunit,*) "WRSI-Parameter: SOS Anomaly file (",&
                trim(temp_filename),") read-in successfully."
!      print *, "WRSI-Parameter: SOS Anomaly file (",&
!                trim(temp_filename),") read-in successfully."
   else
      write(LIS_logunit,*) "ERR MSG:  SOS Anomaly NOT read-in (file: ",&
                           trim(temp_filename),")."
      write(LIS_logunit,*) "End run being called ..."
      call LIS_endrun
   endif

 ! Deallocate arrays:
   call dealloc_arr(real4ptr2dIsNull, real4ptr2d=real4ptr2d)
   call nullify_ptr(real4ptr2dIsNull, real4ptr2d=real4ptr2d)
   deallocate(gCoords)

 end subroutine geowrsi2_readSOS_Bilfile
