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
! !ROUTINE: process_cloud_layers
! \label{process_cloud_layers}
!
! !REVISION HISTORY:
! 09 Jun 2016: James Geiger; Initial specification from pseudo code provided
!                            by Tim Nobis and Mark Conner
! 16 Dec 2021: Eric Kemp: Replaced julhr with YYYYMMDDHH in log.
!
! !INTERFACE:
subroutine process_cloud_layers(cod_layered, base_layered, &
                                top_layered, pcts_layered, &
                                tot_pcts,                  &
                                julhr, cloud_times,        &
                                cod_lmh, pcts_lmh)
!  !USES:
   use LIS_logMod,        only : LIS_logunit

   implicit none
!  !ARGUMENTS:
   integer, parameter :: NC=1024, NR=1024, NL=4

   real, intent(in),  dimension(NL,NC,NR) :: cod_layered,  &
                                             base_layered, &
                                             top_layered,  &
                                             pcts_layered
   real, intent(in), dimension(NC,NR)     :: tot_pcts
   integer, intent(in)                    :: julhr
   real, intent(in),  dimension(NC,NR)    :: cloud_times
   real, intent(out), dimension(3,NC,NR)  :: cod_lmh, pcts_lmh

! !DESCRIPTION:
! This routine processes the four-layer CDFS II cloud optical depth,
! cloud base, cloud top, and cloud percentages into High, Middle, and Low
! cloud optical depths and cloud amounts for computing longwave and shortwave
! radiation.
!
! NOTE: We (Tim Nobis and Mark Conner) are recommending that the
! categorization of a cloud as 'low', 'middle', or 'high' be predicated
! upon the location of the base since the goal of this exercise is
! to characterize impacts of clouds to LW/SW input to the surface.
! The following three if statements will map the four layer data into
! three layers.  The technique uses a weighted
! average to combine optical depths from multiple layers.
!
!  The arguments are:
!  \begin{description}
!  \item[cod\_layered]
!    layered cloud optical depth (one hemisphere)
!  \item[base\_layered]
!    layered cloud base height (one hemisphere)
!  \item[top\_layered]
!    layered cloud top height (one hemisphere)
!  \item[pcts\_layered]
!    layered cloud amount (percentages) (one hemisphere)
!  \item[cod\_lmh]
!    low/middle/high cloud optical depth (one hemisphere)
!  \item[pcts\_lmh]
!    low/middle/high cloud amount (percentages) (one hemisphere)
!  \end{description}
!
!EOP

   integer :: i, j, m, lev
   real    :: weight
   character*10 :: date10_julhr, date10_cloud_time ! EMK

   ! Low    = 1
   ! Middle = 2
   ! High   = 3

   cod_lmh  = 0.0
   pcts_lmh = 0.0

   do j = 1, NR
      do i = 1, NC
         if ( cloud_times(i,j) == 0 ) then
            ! point is not available; skip.
            cycle
         endif
         if ( ( julhr - cloud_times(i,j) .gt. -3 ) .and. &
              ( julhr - cloud_times(i,j) .lt.  6 ) )  then
            do m = 1, NL
               if ( ( pcts_layered(m,i,j) > 0 ) .and. &
                    ( base_layered(m,i,j) > 0 ) ) then

                  if ( base_layered(m,i,j) < 2000 ) then
                     ! identifies any low cloud
                     lev = 1
                  elseif ( base_layered(m,i,j) >= 2000 .and. &
                           base_layered(m,i,j) <= 6000 ) then
                     ! identifies any mid cloud
                     lev = 2
                  elseif ( base_layered(m,i,j) > 6000 ) then
                     ! identifies any high cloud
                     lev = 3
                  endif

                  weight = pcts_lmh(lev,i,j) + pcts_layered(m,i,j)
                  !weight = min(weight, tot_pcts(i,j))
                  if ( weight > 100.0 ) then
                     weight = min(tot_pcts(i,j), 100.0)
                  endif
                  cod_lmh(lev,i,j) = pcts_layered(m,i,j) / weight * &
                                     cod_layered(m,i,j) +           &
                                     pcts_lmh(lev,i,j) / weight *   &
                                     cod_lmh(lev,i,j)
                  pcts_lmh(lev,i,j) = weight
               endif
            enddo
         else
            ! point has invalid time; skip.
            ! EMK Replace julhr output with YYYYMMDDHH
            call AGRMET_julhr_date10 ( julhr, date10_julhr )
            call AGRMET_julhr_date10 ( int(cloud_times(i,j)), &
                 date10_cloud_time )
            write(LIS_logunit,*) '[WARN] cloud_time is not valid'
            write(LIS_logunit,*) '       i,j ', i,j
            !write(LIS_logunit,*) '       julhr', julhr
            write(LIS_logunit,*) '       LIS YYYYMMDDHH ', date10_julhr
            !write(LIS_logunit,*) '       timeLastUpdate', cloud_times(i,j)
            write(LIS_logunit,*) '       timeLastUpdate ', date10_cloud_time
            write(LIS_logunit,*) '       diff ', julhr-cloud_times(i,j)
         endif
      enddo
   enddo

! There was some discussion on the possible existence of a special
! handling path for cloud types identified as 'CB'  (cloud type #1 in CDFS
! speak) and doesn't fill the entire gridcell... this could lead to SW
! reflectance from the cloud side which would reach the ground? If such a
! special path exists, we would recommend using the following logic to
! replace the current CB logic:

! This would be active within an i, j, m loop:

! if ( base_layered(m,i,j) > 0 ) then
!    thickness = top_layered(m,i,j) - base_layered(m,i,j)
!    if ( thickness > 6000 .and. pcts_layered(m,i,j) < 100 ) then
!       treat as if cloud type is CB
!    endif
! endif

end subroutine process_cloud_layers
