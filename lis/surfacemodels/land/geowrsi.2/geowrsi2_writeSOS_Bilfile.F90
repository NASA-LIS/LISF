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
! !ROUTINE:  geowrsi2_writeSOS_Bilfile
!  \label{geowrsi2_writeSOS_Bilfile}
!
! !REVISION HISTORY:
! 31 Jul 2011: Brad Wind; Initial setup
! 15 Jan 2013: KR Arsenault;  Implemented LIS-SOS BIL write routine into LIS
! 25 Oct 2013: KR Arsenault;  Added GeoWRSI2.0 model to LIS-7
!
! !INTERFACE:
 subroutine geowrsi2_writeSOS_Bilfile(n)

! !USES:
  use ESMF
  use LIS_coreMod,    only : LIS_rc, LIS_domain, &
                             LIS_masterproc, LIS_npes
  use LIS_logMod,     only : LIS_logunit, LIS_endrun
  use LIS_fileIOMod,  only : LIS_create_output_directory
  use LIS_historyMod, only : LIS_gather_gridded_output
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use geowrsi2_module
  use geowrsi2_lsmMod
  use geowrsi2_physics_module, only : gMASK_EXC, &
            gWILTING1_DEFAULT, gWILTING2_DEFAULT,& 
            gSOS_INITIAL, gSOS_NOSTART, gSOSCLIM_NA
  use fbil_module

  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n   ! nest
!
! !DESCRIPTION:
! 
!  This routine writes the SOS and SOS anomaly fields into BIL
!  formatted  files.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!EOP

! BIL parameters/variables:
  integer             :: t, c, r, l, s
  integer             :: ftn
  integer             :: index1, ierr
  character(len=LIS_CONST_PATH_LEN) :: dirname
  character(len=4)    :: nest_str
  character(len=4)    :: cyr4
  real*4, allocatable     :: gvar(:)
  real*4, allocatable     :: gtmp(:,:)
  real*4, pointer     :: glb_gvar(:)
  type(charN)         :: outFileFullPathNoExt
  type(FEWSNET_bil__header) :: hdrInf
  type(geowrsi2dec), pointer :: geowrsi2Pt

! ________________________________________________________

   if( LIS_rc%npatch(n,LIS_rc%lsm_index) == 0 ) then
      return
   endif


!- FBIL Write Components:
   if ( LIS_masterproc ) then
      allocate(gCoords)
      allocate(outFileFullPathNoExt%str)
      allocate(glb_gvar(LIS_rc%gnc(n)*LIS_rc%gnr(n)))
   endif
   allocate(gvar(LIS_rc%npatch(n,LIS_rc%lsm_index)))

   if ( LIS_masterproc ) then
      ! Note that the masterproc contains the lower left grid-cell.
      ! If the number of processing elements is 1 then the masterproc
      ! also contains the upper right grid-cell, else the masterproc
      ! must compute the location of the upper right grid-cell.
      !
      ! Note that the BIL library expects grid-cells to be defined
      ! by their upper left corner.  LIS defines grid-cells on their
      ! centers.  When setting gCoords, nudge the values to conform
      ! to the BIL library.
      if ( LIS_npes == 1 ) then
         gCoords%minLat = LIS_rc%gridDesc(n,4) + LIS_rc%gridDesc(n,10)/2.0
         gCoords%maxLat = LIS_rc%gridDesc(n,7) + LIS_rc%gridDesc(n,10)/2.0
         gCoords%minLon = LIS_rc%gridDesc(n,5) - LIS_rc%gridDesc(n,9)/2.0
         gCoords%maxLon = LIS_rc%gridDesc(n,8) - LIS_rc%gridDesc(n,9)/2.0
         gCoords%pixLat = LIS_rc%gridDesc(n,10)
         gCoords%pixLon = LIS_rc%gridDesc(n,9)
      else
         gCoords%minLat = LIS_rc%gridDesc(n,4) + LIS_rc%gridDesc(n,10)/2.0
         gCoords%maxLat = ( LIS_rc%gridDesc(n,4) +                      &
                            (LIS_rc%gnr(n)-1)*LIS_rc%gridDesc(n,10) ) + &
                          LIS_rc%gridDesc(n,10)/2.0
         gCoords%minLon = LIS_rc%gridDesc(n,5) - LIS_rc%gridDesc(n,9)/2.0
         gCoords%maxLon = ( LIS_rc%gridDesc(n,5) +                     &
                            (LIS_rc%gnc(n)-1)*LIS_rc%gridDesc(n,9) ) - &
                          LIS_rc%gridDesc(n,9)/2.0
         gCoords%pixLat = LIS_rc%gridDesc(n,10)
         gCoords%pixLon = LIS_rc%gridDesc(n,9)
      endif

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

      write(hdrInf%byteorder, '(a)') 'I' ! default value
      write(hdrInf%layout, '(a)') 'BIL'  ! default value
     ! hdrInf%nbits     = ! deferred until write time
      hdrInf%xdim      = gCoords%pixLon
      hdrInf%ydim      = gCoords%pixLat
      hdrInf%ncols     = gCoords%maxX - gCoords%minX + 1
      hdrInf%nrows     = gCoords%maxY - gCoords%minY + 1
      hdrInf%nbands    = 1               ! default value
      hdrInf%ulxmap    = gCoords%minLon
      hdrInf%ulymap    = gCoords%maxLat
   endif

 ! Growing season count:
   s = geowrsi2_struc(n)%season_count

 ! Write name of nest:
   write(unit=nest_str, fmt='(a2,i2.2)') '.d',n

 ! =====  SOS  ======

   if ( LIS_masterproc ) then
      write(unit=cyr4,fmt='(i4.4)') LIS_rc%yr
      dirname = trim(LIS_rc%odir)
      write(outFileFullPathNoExt%str, '(a)') &
               trim(dirname)//"/"//"SOS_current_"//cyr4//nest_str
   endif

   gvar = LIS_rc%udef
   do t = 1, LIS_rc%npatch(n,LIS_rc%lsm_index)
      gvar(t) = geowrsi2_struc(n)%wrsi(t)%sos_write(s) ! CHANGE TO ACCOMMODATE
                                                       ! FOR MULTI-GROWING SEASONS
   enddo

   call LIS_gather_gridded_output(n, LIS_rc%lsm_index, gtmp, gvar)

   if ( LIS_masterproc ) then
    ! Write data out:
      glb_gvar = reshape(gtmp, (/LIS_rc%gnc(n)*LIS_rc%gnr(n)/))
      hdrInf%nbits  = gBIL_INT4;
      call write_FEWSNET_bil(outFileFullPathNoExt, hdrInf,&
                             real4ptr1d=glb_gvar)
      deallocate(gtmp)
   endif

 ! =====  SOSa  ======

   if ( LIS_masterproc ) then
      dirname = trim(LIS_rc%odir)
      write(outFileFullPathNoExt%str, '(a)') &
               trim(dirname)//"/"//"SOS_anomaly_"//cyr4//nest_str
   endif

   gvar = LIS_rc%udef
   do t = 1, LIS_rc%npatch(n,LIS_rc%lsm_index)
      gvar(t) = geowrsi2_struc(n)%wrsi(t)%sosa_write(s) ! CHANGE TO ACCOMMODATE
                                                        ! FOR MULTI-GROWING SEASONS
   enddo

   call LIS_gather_gridded_output(n, LIS_rc%lsm_index, gtmp, gvar)

   if ( LIS_masterproc ) then
    ! Write data out:
      glb_gvar = reshape(gtmp, (/LIS_rc%gnc(n)*LIS_rc%gnr(n)/))
      hdrInf%nbits  = gBIL_INT4;
      call write_FEWSNET_bil(outFileFullPathNoExt, hdrInf,&
                             real4ptr1d=glb_gvar)
      deallocate(gtmp)
   endif

 ! Reinitialize SOS/SOSa fields now written:
!   print *, " Reinitializing SOS/SOSa after writing SOS/SOSa files ..."
   do t = 1, LIS_rc%npatch(n,LIS_rc%lsm_index)
      geowrsi2Pt => geowrsi2_struc(n)%wrsi(t)
    ! Calculate SOS only over unmasked areas, omitting areas where there is a mask.
      if( (geowrsi2Pt%Mask == CInt2(real(geowrsi2_udef, 8))) .or. &
          (geowrsi2Pt%Mask == gMASK_EXC) ) then   ! 0
         geowrsi2Pt%SOS = gSOS_INITIAL   ! 0
      else
         geowrsi2Pt%SOS = gSOS_NOSTART   ! 60
      endif
      geowrsi2Pt%SOSa     = gSOSCLIM_NA  ! 0
      geowrsi2Pt%Wilting1 = gWILTING1_DEFAULT
      geowrsi2Pt%Wilting2 = gWILTING2_DEFAULT
   end do


   if ( LIS_masterproc ) then
      deallocate(gCoords)
      deallocate(outFileFullPathNoExt%str)
      deallocate(glb_gvar)
   endif
   deallocate(gvar)

 end subroutine geowrsi2_writeSOS_Bilfile

