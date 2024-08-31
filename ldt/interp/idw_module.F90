!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! by Matthew Garcia, GEST Research Associate
!    Hydrological Sciences Branch (Code 614.3)
!    NASA-GSFC 
!
! Calculation of inverse-distance weights
! -- Cartesian distance calculation      -- functions cartdist2,3
! -- inverse-distance weight calculation -- functions idweight2,3
! -- quadrant-based IDW                  -- subroutine idwquads
! -- ranked nearest-neighbors IDW        -- subroutine idwranked
! -- all-neighbors IDW                   -- subroutine idwalln
! -- normalization of gridded weights    -- subroutine normalize
!
! Version History
!   12Aug04  Matthew Garcia  Program idwgrid (StIM v0.9) delivered to ARMS project
!   21Aug04  Matthew Garcia  IDW calculation processes distributed to invdist.F90
!   20Apr05  Matthew Garcia  Added functions for 3D interpolation
!
module idw_module
!
  implicit none
!
contains
!
!
real(4) function cartdist2(stnx,stny,x,y)
!
! argument variables
  real(4), intent(IN) :: stnx,stny,x,y
!
! local variables
  real(4) :: xdist,ydist,dsq
!
  xdist = stnx - x
  ydist = stny - y
  dsq = (xdist * xdist) + (ydist * ydist)
  cartdist2 = sqrt(dsq)
end function cartdist2
!
!
real(4) function cartdist3(stnx,stny,stnz,x,y,z)
!
! argument variables
  real(4), intent(IN) :: stnx,stny,stnz,x,y,z
!
! local variables
  real(4) :: xdist,ydist,zdist,dsq
!
  xdist = stnx - x
  ydist = stny - y
  zdist = stnz - z
  dsq = (xdist * xdist) + (ydist * ydist) + (zdist * zdist)
  cartdist3 = sqrt(dsq)
end function cartdist3
!
!
real(4) function idweight2(stnx,stny,x,y,ord)
!
! argument variables
  real(4), intent(INOUT) :: stnx,stny,x,y
  integer, intent(IN) :: ord
!
! local variables
  real(4) :: dist
!
  dist = cartdist2(stnx,stny,x,y)
  if (dist.eq.0) then 
        idweight2 = 1.0 
  else 
        idweight2 = dist ** (-1 * ord)
  endif
end function idweight2
!
!
real(4) function idweight3(stnx,stny,stnz,x,y,z,ord)
!
! argument variables
  real(4), intent(INOUT) :: stnx,stny,stnz,x,y,z
  integer, intent(IN) :: ord
!
! local variables
  real(4) :: dist
!
  dist = cartdist3(stnx,stny,stnz,x,y,z)
  if (dist.eq.0) then 
        idweight3 = 1.0 
  else 
        idweight3 = dist ** (-1 * ord)
  endif
end function idweight3
!
!
subroutine idwquads(ns,nn,nc,nr,lx,ry,inc,ord,sdata,stnarr,wgtarr)
!
! argument variables
  integer, intent(IN) :: ns,nn,nc,nr
  real(4), intent(IN) :: lx,ry,inc
  integer, intent(INOUT) :: ord
  real(4), intent(INOUT) :: sdata(:,:)
  integer(4), intent(INOUT) :: stnarr(:,:,:) ! station numbers for each grid cell
  real(4), intent(INOUT) :: wgtarr(:,:,:) ! station weights for each grid cell 
!
! local variables
  real(4) :: xloc,yloc
  real(4) :: weight
  integer :: s,i,j
!
  print *,'MSG: idw_module -- establishing quadrant-based non-normal weights.'
  do s = 1,ns
    do i = 1,nc
      xloc = lx + inc * (i - 1)
          do j = 1,nr
            yloc = ry - inc * (j - 1)
                weight = idweight2(sdata(2,s),sdata(3,s),xloc,yloc,ord)
            if (weight.eq.1) then ! spot-on
              stnarr(nn+1,i,j) = sdata(1,s)
            else
              if (sdata(2,s).ge.xloc.and.sdata(3,s).gt.yloc) then ! NE (Q1)
                if (weight.gt.wgtarr(1,i,j)) then
                  wgtarr(1,i,j) = weight
                  stnarr(1,i,j) = sdata(1,s)
                endif
              else if (sdata(2,s).gt.xloc.and.sdata(3,s).le.yloc) then ! SE (Q2)
                if (weight.gt.wgtarr(2,i,j)) then
                  wgtarr(2,i,j) = weight
                  stnarr(2,i,j) = sdata(1,s)
                endif
              else if (sdata(2,s).le.xloc.and.sdata(3,s).lt.yloc) then ! SW (Q3)
                if (weight.gt.wgtarr(3,i,j)) then
                  wgtarr(3,i,j) = weight
                  stnarr(3,i,j) = sdata(1,s)
                endif
              else if (sdata(2,s).lt.xloc.and.sdata(3,s).ge.yloc) then ! NW (Q4)
                if (weight.gt.wgtarr(4,i,j)) then
                  wgtarr(4,i,j) = weight
                  stnarr(4,i,j) = sdata(1,s)
                endif
              endif
            endif
      end do
    end do
  end do
!
end subroutine idwquads
!
!
subroutine idwranked(ns,nn,nc,nr,lx,ry,inc,ord,sdata,stnarr,wgtarr)
!
! argument variables
  integer, intent(IN) :: ns,nn,nc,nr
  real(4), intent(IN) :: lx,ry,inc
  integer, intent(INOUT) :: ord
  real(4), intent(INOUT) :: sdata(:,:)
  integer(4), intent(INOUT) :: stnarr(:,:,:) ! station numbers for each grid cell
  real(4), intent(INOUT) :: wgtarr(:,:,:) ! station weights for each grid cell 
!
! local variables
  real(4) :: xloc,yloc
  real(4) :: weight
  integer :: n,s,i,j
!
  print *,'MSG: idw_module -- establishing ranked non-normal weights.'
  do s = 1,ns
    do i = 1,nc
      xloc = lx + inc * (i - 1)
      do j = 1,nr
            yloc = ry - inc * (j - 1)
            weight = idweight2(sdata(2,s),sdata(3,s),xloc,yloc,ord)
            if (weight.eq.1) then ! spot-on
                  stnarr(nn+1,i,j) = sdata(1,s)
                  wgtarr(1,i,j) = weight
                  stnarr(1,i,j) = sdata(1,s)
                 do n = 2,nn
                   wgtarr(n,i,j) = 0.0
                   stnarr(n,i,j) = 0
                 end do
            else
                 if (wgtarr(1,i,j).ne.1) then
                   do n = nn-1,1,-1
                     if (weight.gt.wgtarr(n,i,j)) then
                       wgtarr(n+1,i,j) = wgtarr(n,i,j)
                       stnarr(n+1,i,j) = stnarr(n,i,j)
                       wgtarr(n,i,j) = weight
                       stnarr(n,i,j) = sdata(1,s)
                     end if
                   end do
            end if
        end if
      end do
    end do
  end do
!
end subroutine idwranked
!
!
subroutine idwalln(ns,nc,nr,lx,ry,inc,ord,sdata,wgtarr)
!
! argument variables
  integer, intent(IN) :: ns,nc,nr
  real(4), intent(IN) :: lx,ry,inc
  integer, intent(INOUT) :: ord
  real(4), intent(INOUT) :: sdata(:,:)
  real(4), intent(INOUT) :: wgtarr(:,:,:) ! station weights for each grid cell 
!
! local variables
  real(4) :: xloc,yloc
  real(4) :: weight,maxwgt
  integer :: n,s,i,j
!
  print *,'MSG: idw_module -- establishing all-stations non-normal weights.'
  do s = 1,ns
    do i = 1,nc
      xloc = lx + inc * (i - 1)
        do j = 1,nr
            yloc = ry - inc * (j - 1)
            weight = idweight2(sdata(2,s),sdata(3,s),xloc,yloc,ord)
            if (weight.eq.1) then ! spot-on
              do n = 1,ns
                wgtarr(n,i,j) = 0.0
              end do
              wgtarr(s,i,j) = weight
            else
              maxwgt = 0.0
              do n = 1,ns
                if (wgtarr(n,i,j).gt.maxwgt) maxwgt = wgtarr(n,i,j)
              end do          
              if (maxwgt.ne.1) then
                wgtarr(s,i,j) = idweight2(sdata(2,s),sdata(3,s),xloc,yloc,ord)
              end if
            end if
      end do
    end do
  end do
!
end subroutine idwalln
!
!
subroutine normalize(nn,nc,nr,wgtarr)
!
! argument variables
  integer, intent(IN) :: nn,nc,nr
  real(4), intent(INOUT) :: wgtarr(:,:,:) ! station weights for each grid cell 
!
! local variables
  integer :: s,i,j
!
  print *,'MSG: idw_module -- establishing gridded normalized weights for stations.'
  do s = 1,nn
    do i = 1,nc
      do j = 1,nr
        wgtarr(nn+1,i,j) = wgtarr(nn+1,i,j) + wgtarr(s,i,j)
      end do
    end do
  end do
  do s = 1,nn
    do i = 1,nc
      do j = 1,nr
        wgtarr(s,i,j) = wgtarr(s,i,j) / wgtarr(nn+1,i,j) 
      end do
    end do
  end do
!
end subroutine normalize
!
!
end module idw_module
