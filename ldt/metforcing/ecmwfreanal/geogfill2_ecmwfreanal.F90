!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! 
! !ROUTINE: geogfill2_ecmwfreanal
!  \label{geogfill2_ecmwfreanal}
!
! !REVISION HISTORY:
!  20 Sep 2002: Urszula Jambor; Modified original geogfill routine
!               for use with Reanalysis forcing sets prepared by
!               Aaron Berg via NSIPP
!  29 Jan 2003: Urszula Jambor; Removed LDAS & GRID modules from list of
!               arguments, instead pass only needed variables directly
!               (nc,nr,fimask).
!  27 Jan 2004: Matt Rodell; Make sure all points are filled, rename
!               data array to fdata.
! !INTERFACE:
subroutine geogfill2_ecmwfreanal(n, nc,nr,geogmask,fdata,v,tmask)
! !USES:
  use LDT_coreMod, only : LDT_domain

  implicit none
! !ARGUMENTS:   
  integer, intent(in) :: n 
  integer             :: nc
  integer             :: nr
  logical*1           :: geogmask(nc,nr)    
  real                :: fdata(nc,nr)       
  integer             :: v
  logical*1           :: tmask(nc,nr)
!
! !DESCRIPTION:
!  Fill in grid points that are invalid in LDT due to differences in
!  geography between forcing and land surface.
!
!  Based on original geogfill.F90 for use with Reanalysis forcing 
!  fdata, to account for differences in land masks.  
!  NOTE** for Reanalysis fdata, v=6 is the V-component of wind,
!  which is always set to ZERO since the U-component (v=5) is assigned 
!  the absolute value of wind speed.
!
!  For v=1:temperature, v=2:specific humidity, v=4:LW radiation,
!  v=5:wind, and v=7:pressure, data values of zero are not allowed to 
!  contribute to average fill-value.
!  For v=3:SW radiation, v=8,9:precipitation, use the mask defined by 
!  valid temperature data points to establish whether to seek fill-value
!  for given column and row, rather than include vs exclude zeroes.
!
!  The arguments are:  
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[nc]
!   number of columns of the grid (along the east-west dimension)
!  \item[nr]
!   number of rows of the grid (along the north-south  dimension)
!  \item[geogmask]
!   interpolated forcing bitmap
!  \item[fdata]
!   interpolated forcing
!  \item[v]
!   index of the variable
!  \item[tmask]
!   valid temperature points (used as the land mask)
!  \end{description}
!EOP
  integer :: i,j,kk,kt,ii,jj          ! loop counters

  integer :: smask(nc,nr)                ! grid points where ldas
  real    :: ssum,savg                      ! sum, average--fill in points
  integer :: scnt                        ! counter--fill in points
  logical*1 :: filled 

  !=== initializing sub-surface parameter bitmap ===========================
  smask = 0

  !=== defining sub-surface parameter bitmap ===============================
  do j=1,nr
     do i=1,nc
        if (LDT_domain(n)%gindex(i,j).ne. -1 .and. &
             .not.(geogmask(i,j))) then
           smask(i,j) = 1
        endif
     enddo
  enddo

  !=== filling in points where sub-surface parameter bitmap is equal to 1
  do j = 1, nr
     do i = 1, nc
        ssum = 0
        scnt = 0
        if (smask(i,j) .gt. 0) then              ! problem land point
           if (j .le. 2 .or. j .ge. nr-1) then
              fdata(i,j) = fdata(i,j)
           elseif (i .le. 2) then                 ! first two columns
              if(i .eq. 1) then
                 kt = nc-2
              else
                 kt = nc-1
              endif
              do kk = 1, 5
                 kt = kt + 1
                 do jj = j-2, j+2
                    if ( (v.ne.3) .and. (v.lt.8) ) then
                       if ((fdata(kt,jj).gt.0.0) .and. fdata(kt,jj).ne.-1) then
                          ssum = ssum + fdata(kt,jj)
                          scnt = scnt + 1
                       endif
                    else if ( tmask(i,j) .and. fdata(kt,jj).ne.-1) then
                       if (smask(kt,jj).eq.0) then
                          ssum = ssum + fdata(kt,jj)
                          scnt = scnt + 1
                       endif
                    endif
                 enddo
                 if (kt .eq. nc) kt = 0
              enddo
              if (scnt .ne. 0) then
                 savg = ssum / scnt      ! average of surrounding points
              else
                 savg = -1               ! set flag for later
              endif
              fdata(i,j) = savg                 ! assign fill-in fdata point
           elseif (i .ge. nc-1) then     ! last two columns
              kt = i-3
              do kk = 1, 5
                 kt = kt + 1
                 do jj = j-2, j+2
                    if ( (v.ne.3) .and. (v.lt.8) ) then
                       if ((fdata(kt,jj).gt.0.0) .and. fdata(kt,jj).ne.-1) then
                          ssum = ssum + fdata(kt,jj)
                          scnt = scnt + 1
                       endif
                    else if ( tmask(i,j) .and. fdata(kt,jj).ne.-1 ) then
                       if (smask(kt,jj).eq.0) then
                          ssum = ssum + fdata(kt,jj)
                          scnt = scnt + 1
                       endif
                    endif
                 enddo
                 if (kt .eq. nc) kt = 0
              enddo
              if (scnt .ne. 0) then
                 savg = ssum / scnt
              else
                 savg = -1
              endif
              fdata(i,j) = savg
           else                                   ! all other points
              do ii = i-5, i+5
                 do jj = j-5, j+5
                    if ( (v.ne.3) .and. (v.lt.8) ) then
                       if ((fdata(ii,jj).gt.0.0) .and. fdata(ii,jj).ne.-1 &
                            .and.ii.gt.0.and.jj.gt.0) then
                          ssum = ssum + fdata(ii,jj)
                          scnt = scnt + 1
                       endif
                    else if ( tmask(i,j) .and. fdata(ii,jj).ne.-1 ) then
                       if (smask(ii,jj).eq.0 .and.ii.gt.0.and.jj.gt.0) then
                          ssum = ssum + fdata(ii,jj)
                          scnt = scnt + 1
                       endif
                    endif
                 enddo
              enddo
              if (scnt .ne. 0) then
                 savg = ssum / scnt
              else
                 savg = -1
              endif
              fdata(i,j) = savg
           endif
        endif
     enddo
  enddo

  if (v==1) then ! assign tmask for future use
     do j=1,nr
        do i=1,nc
           if (fdata(i,j) > 0.0) then
              tmask(i,j) = .true.
           endif
        enddo
     enddo
  endif

  !     fill in remaining points with the nearest valid forcing value in 
  !      the direction of the most data (largest grid).
  do j=1,nr
     do i=1,nc
        if (smask(i,j) .eq. 1) then
           if (((fdata(i,j) .lt. -0.99).and.(fdata(i,j) .gt. -1.01)) & 
                .or. ((j .le. 2) .or. (j .ge. (nr-1)))) then
              filled = .false.
              if (j .le. nr/2) then
                 do jj = j,nr
                    if (i .le. nc/2) then
                       do ii=i,nc
                          if (tmask(ii,jj)) then
                             fdata(i,j) = fdata(ii,jj)
                             filled = .true.
                             exit !out of ii do loop
                          end if
                       end do
                    else
                       do ii=i,1,-1
                          if (tmask(ii,jj)) then
                             fdata(i,j) = fdata(ii,jj)
                             filled = .true.
                             exit !out of ii do loop
                          end if
                       end do
                    end if
                    if (filled) exit !out of jj do loop
                 end do
              else
                 do jj = j,1,-1
                    if (i .le. nc/2) then
                       do ii=i,nc
                          if (tmask(ii,jj))then
                             fdata(i,j) = fdata(ii,jj)
                             filled = .true.
                             exit !out of ii do loop
                          end if
                       end do
                    else
                       do ii=i,1,-1
                          if (tmask(ii,jj)) then
                             fdata(i,j) = fdata(ii,jj)
                             filled = .true.
                             exit !out of ii do loop
                          end if
                       end do
                    end if
                    if (filled) exit !out of jj do loop
                 end do
              end if !j
           end if !fdata
        end if !smask
     end do !i
  end do !j

  return
end subroutine geogfill2_ecmwfreanal

