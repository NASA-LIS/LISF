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
! !ROUTINE: interp_agrmetvar
!  \label{interp_agrmetvar}
! 
! !INTERFACE:
subroutine interp_agrmetvar(n,ip,gi,udef,varfield,imax,jmax)
!
! 31 MAR 2010 a bunch of modifications for handling various combinations of 
!             different resolution datasets..........Michael Shaw/WXE
!
! !USES:
  use LIS_coreMod,       only : LIS_rc
  use AGRMET_forcingMod, only : agrmet_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n,imax,jmax 
  integer, intent(in)    :: ip
  real, intent(in)       :: gi(2,imax,jmax)
  real, intent(in)       :: udef
  real, intent(inout)    :: varfield(LIS_rc%lnc(n),LIS_rc%lnr(n))
!
! !DESCRIPTION:
!   This subroutine interpolates a given AFWA field 
!   to the LIS domain (from polar stereographic grid to the LIS grid)
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[ip]
!    interpolation option
!  \item[gi]
!    input AGRMET field (both hemispheres)
!  \item[udef]
!    undefined value in the input field
!  \item[varfield]
!    interpolated field in the LIS grid
!  \end{description}
!
!  The routines invoked are: 
!
!  \begin{description}
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[neighbor\_interp](\ref{neighbor_interp}) \newline
!    spatially interpolate the forcing data using neighbor interpolation
! \end{description}
!EOP
  real,allocatable                   :: gi_temp(:)
  integer                :: end
  logical*1,allocatable              :: li(:)
  integer                :: ihemi,hemi,i,j,iret,totdim
  real                   :: gridDesco(200)
  logical*1, allocatable :: lo_nh(:),lo_sh(:)
  real, allocatable      :: go_nh(:),go_sh(:)
  real, allocatable      :: comb_data(:)
  real, allocatable      :: varfield_temp(:,:),shfield_temp(:,:)
  integer                :: c,r,o,p,jj,jadj,jb,je,jindex
!  real                   :: nhfield(agrmet_struc(n)%hemi_nc(1),agrmet_struc(n)%hemi_nr(1))
!  real                   :: shfield(agrmet_struc(n)%hemi_nc(2),agrmet_struc(n)%hemi_nr(2))
  real, allocatable      :: nfield(:,:)
 !

   if(agrmet_struc(n)%interp.gt.0) then
!<kluge -- cod testing>
!      if(ip.eq.1) then 
      if(ip.eq.1.or.ip.eq.3) then 
!</kluge -- cod testing>
         if(agrmet_struc(n)%gridspan.eq.1) then 
            allocate(lo_nh(agrmet_struc(n)%mo1))
            allocate(go_nh(agrmet_struc(n)%mo1))
         elseif(agrmet_struc(n)%gridspan.eq.2) then 
            allocate(lo_sh(agrmet_struc(n)%mo2))
            allocate(go_sh(agrmet_struc(n)%mo2))
         else
            allocate(lo_nh(agrmet_struc(n)%mo1))
            allocate(lo_sh(agrmet_struc(n)%mo2))
            allocate(go_nh(agrmet_struc(n)%mo1))
            allocate(go_sh(agrmet_struc(n)%mo2))
         endif
         alloCate(comb_data(agrmet_struc(n)%mo1+agrmet_struc(n)%mo2))
         do ihemi = agrmet_struc(n)%shemi,agrmet_struc(n)%nhemi
            if(agrmet_struc(n)%global_or_hemi==1)then !only do this for global data
               jb=int(agrmet_struc(n)%global_or_hemi*(2-ihemi)*int(jmax/2)+1)
               je=int(int(jmax)+agrmet_struc(n)%global_or_hemi*(1-ihemi)*int(jmax/2))
               jadj=agrmet_struc(n)%global_or_hemi*int(jmax/2)*(2-ihemi)
            else
               jb=1
               je=jmax
               jadj=0
            endif 

! for nonglobal data, use the "2" structures if this is an interp for nonstandard input
! Also, li and gi_temp should be allocatable since they'll be different depending on the
! size of any nonstandard input (global or just different from basis agrmet grid)

            if(agrmet_struc(n)%global_or_hemi == 0)then
               hemi=ihemi
               allocate(gi_temp(imax*jmax),li(imax*jmax))
               totdim=imax*jmax
            else
 
! for global input, going to split the data into hemispheres for most direct comparison 
! to the old polar stereo grids (in Totalview, e.g.) and will always be using the "2" structures 
! since this is definitely different from "Standard" polar stereo native grid

               hemi=1
               allocate(gi_temp(int(imax*jmax/2)),li(int(imax*jmax/2)))
               totdim=int(imax*jmax/2)
            endif
            gi_temp=LIS_rc%udef
            li = .false.


            do i=1,imax
               do j=jb,je
                  if(gi(hemi,i,j).ne.udef) then 
                     jj = j-jadj
                     li(i+(jj-1)*imax) = .true.
                     gi_temp(i+(jj-1)*imax) = gi(hemi,i,j)
                  endif
               enddo
            enddo

! just make a field to check the array if diff grids (in Totalview, e.g.) - get rid of it immediately - don't need it for anything else

            if(agrmet_struc(n)%diff_grid == 1)then
               if(agrmet_struc(n)%global_or_hemi == 1)then 
                  allocate(nfield(imax,int(jmax/2)))
                  nfield=RESHAPE(gi_temp(1:int(imax*jmax/2)),(/imax,int(jmax/2)/))
               else
                  allocate(nfield(imax,jmax))
                  nfield=RESHAPE(gi_temp(1:int(imax*jmax)),(/imax,jmax/))
               endif
               deallocate(nfield)
            endif

            gridDesco(1) = LIS_rc%gridDesc(n,1) 
            gridDesco(2) = agrmet_struc(n)%hemi_nc(ihemi)
            gridDesco(3) = agrmet_struc(n)%hemi_nr(ihemi)
            gridDesco(5) = LIS_rc%gridDesc(n,5)
            gridDesco(8) = LIS_rc%gridDesc(n,8)
            gridDesco(6) = LIS_rc%gridDesc(n,6)
            gridDesco(9) = LIS_rc%gridDesc(n,9)
            gridDesco(10) = LIS_rc%gridDesc(n,10)
            if(LIS_rc%gridDesc(n,1).ne.0)gridDesco(11) = LIS_rc%gridDesc(n,11)
            !YDT 10/2/07
            !gridDesco(20) = 255
            gridDesco(20) = 0 
            if(agrmet_struc(n)%gridspan.eq.1.or.agrmet_struc(n)%gridspan.eq.2) then 
               gridDesco(4) = LIS_rc%gridDesc(n,4)
               gridDesco(7) = LIS_rc%gridDesc(n,7)
            else
               if(ihemi.eq.1) then 
                  gridDesco(4) = LIS_rc%gridDesc(n,9)/2
                  gridDesco(7) = LIS_rc%gridDesc(n,7)
               else
                  gridDesco(4) = LIS_rc%gridDesc(n,4)
                  gridDesco(7) = -LIS_rc%gridDesc(n,9)/2
               endif
            endif
            if(ihemi.eq.1) then

! Need to do a different interp for nonstandard input with "2" structures

!<kluge -- cod testing>
               if ( ip == 3 ) then
                  call bilinear_interp(gridDesco,li,gi_temp,&
                       lo_nh,go_nh,totdim,agrmet_struc(n)%mo1,&
                       agrmet_struc(n)%rlat1_nh4,agrmet_struc(n)%rlon1_nh4,&
                       agrmet_struc(n)%w111_nh4,agrmet_struc(n)%w121_nh4,&
                       agrmet_struc(n)%w211_nh4,agrmet_struc(n)%w221_nh4,&
                       agrmet_struc(n)%n111_nh4,agrmet_struc(n)%n121_nh4,&
                       agrmet_struc(n)%n211_nh4,agrmet_struc(n)%n221_nh4,LIS_rc%udef,&
                       iret)
               elseif(imax == agrmet_struc(n)%imax3)then
!               if(imax == agrmet_struc(n)%imax3)then
!</kluge -- cod testing>
                  call bilinear_interp(gridDesco,li,gi_temp,&
                       lo_nh,go_nh,totdim,agrmet_struc(n)%mo1,&
                       agrmet_struc(n)%rlat1_nh3,agrmet_struc(n)%rlon1_nh3,&
                       agrmet_struc(n)%w111_nh3,agrmet_struc(n)%w121_nh3,&
                       agrmet_struc(n)%w211_nh3,agrmet_struc(n)%w221_nh3,&
                       agrmet_struc(n)%n111_nh3,agrmet_struc(n)%n121_nh3,&
                       agrmet_struc(n)%n211_nh3,agrmet_struc(n)%n221_nh3,LIS_rc%udef,&
                       iret)
               elseif(imax == agrmet_struc(n)%imax2)then 
                  call bilinear_interp(gridDesco,li,gi_temp,&
                       lo_nh,go_nh,totdim,agrmet_struc(n)%mo1,&
                       agrmet_struc(n)%rlat1_nh2,agrmet_struc(n)%rlon1_nh2,&
                       agrmet_struc(n)%w111_nh2,agrmet_struc(n)%w121_nh2,&
                       agrmet_struc(n)%w211_nh2,agrmet_struc(n)%w221_nh2,&
                       agrmet_struc(n)%n111_nh2,agrmet_struc(n)%n121_nh2,&
                       agrmet_struc(n)%n211_nh2,agrmet_struc(n)%n221_nh2,LIS_rc%udef,&
                       iret)
               else
                  call bilinear_interp(gridDesco,li,gi_temp,&
                       lo_nh,go_nh,totdim,agrmet_struc(n)%mo1,&
                       agrmet_struc(n)%rlat1_nh,agrmet_struc(n)%rlon1_nh,&
                       agrmet_struc(n)%w111_nh,agrmet_struc(n)%w121_nh,&
                       agrmet_struc(n)%w211_nh,agrmet_struc(n)%w221_nh,&
                       agrmet_struc(n)%n111_nh,agrmet_struc(n)%n121_nh,&
                       agrmet_struc(n)%n211_nh,agrmet_struc(n)%n221_nh,LIS_rc%udef,&
                       iret)
               endif
            elseif(ihemi.eq.2) then
!<kluge -- cod testing>
               if ( ip == 3 ) then
                  call bilinear_interp(gridDesco,li,gi_temp,&
                       lo_sh,go_sh,totdim,agrmet_struc(n)%mo2,&
                       agrmet_struc(n)%rlat1_sh4,agrmet_struc(n)%rlon1_sh4,&
                       agrmet_struc(n)%w111_sh4,agrmet_struc(n)%w121_sh4,&
                       agrmet_struc(n)%w211_sh4,agrmet_struc(n)%w221_sh4,&
                       agrmet_struc(n)%n111_sh4,agrmet_struc(n)%n121_sh4,&
                       agrmet_struc(n)%n211_sh4,agrmet_struc(n)%n221_sh4,LIS_rc%udef,&
                       iret)
               elseif(imax == agrmet_struc(n)%imax3)then
!               if(imax == agrmet_struc(n)%imax3)then 
!</kluge -- cod testing>
                  call bilinear_interp(gridDesco,li,gi_temp,&
                       lo_sh,go_sh,totdim,agrmet_struc(n)%mo2,&
                       agrmet_struc(n)%rlat1_sh3, agrmet_struc(n)%rlon1_sh3,&
                       agrmet_struc(n)%w111_sh3,agrmet_struc(n)%w121_sh3,&
                       agrmet_struc(n)%w211_sh3,agrmet_struc(n)%w221_sh3,&
                       agrmet_struc(n)%n111_sh3,agrmet_struc(n)%n121_sh3,&
                       agrmet_struc(n)%n211_sh3,agrmet_struc(n)%n221_sh3,&
                       LIS_rc%udef,iret)
               elseif(imax == agrmet_struc(n)%imax2)then
                  call bilinear_interp(gridDesco,li,gi_temp,&
                       lo_sh,go_sh,totdim,agrmet_struc(n)%mo2,&
                       agrmet_struc(n)%rlat1_sh2, agrmet_struc(n)%rlon1_sh2,&
                       agrmet_struc(n)%w111_sh2,agrmet_struc(n)%w121_sh2,&
                       agrmet_struc(n)%w211_sh2,agrmet_struc(n)%w221_sh2,&
                       agrmet_struc(n)%n111_sh2,agrmet_struc(n)%n121_sh2,&
                       agrmet_struc(n)%n211_sh2,agrmet_struc(n)%n221_sh2,&
                       LIS_rc%udef,iret)
               else
                  call bilinear_interp(gridDesco,li,gi_temp,&
                       lo_sh,go_sh,totdim,agrmet_struc(n)%mo2,&
                       agrmet_struc(n)%rlat1_sh, agrmet_struc(n)%rlon1_sh,&
                       agrmet_struc(n)%w111_sh,agrmet_struc(n)%w121_sh,&
                       agrmet_struc(n)%w211_sh,agrmet_struc(n)%w221_sh,&
                       agrmet_struc(n)%n111_sh,agrmet_struc(n)%n121_sh,&
                       agrmet_struc(n)%n211_sh,agrmet_struc(n)%n221_sh,&
                       LIS_rc%udef,iret)
               endif
            endif
            deallocate(li,gi_temp)
         enddo
         if(agrmet_struc(n)%gridspan.eq.1) then
            comb_data(1:agrmet_struc(n)%mo2) = LIS_rc%udef
            comb_data(agrmet_struc(n)%mo2+1:agrmet_struc(n)%mo1+agrmet_struc(n)%mo2) = go_nh(:)
         elseif(agrmet_struc(n)%gridspan.eq.2) then
            comb_data(1:agrmet_struc(n)%mo2) = go_sh(:)
            comb_data(agrmet_struc(n)%mo2+1:agrmet_struc(n)%mo1+agrmet_struc(n)%mo2) = LIS_rc%udef
         else
            comb_data(1:agrmet_struc(n)%mo2) = go_sh(:)
            comb_data(agrmet_struc(n)%mo2+1:agrmet_struc(n)%mo1+agrmet_struc(n)%mo2) = go_nh(:)
         endif
! nhfield and shfield simply for debugging (view of arrays) - low overhead so keep
! SVK - The following lines were introduced by AFWA and are problematic when running in 
! parallel. They are turned off as they are purely for debugging purposes. 
!         nhfield = RESHAPE(go_nh(1:agrmet_struc(n)%mo1),(/agrmet_struc(n)%hemi_nc(1),agrmet_struc(n)%hemi_nr(1)/))
!         shfield = RESHAPE(go_sh(1:agrmet_struc(n)%mo2),(/agrmet_struc(n)%hemi_nc(2),agrmet_struc(n)%hemi_nr(2)/))
         varfield = RESHAPE(comb_data(1:agrmet_struc(n)%mo1+agrmet_struc(n)%mo2),(/LIS_rc%lnc(n),LIS_rc%lnr(n)/))
         if(agrmet_struc(n)%gridspan.eq.1) then 
            deallocate(lo_nh)
            deallocate(go_nh)
         elseif(agrmet_struc(n)%gridspan.eq.2) then 
            deallocate(lo_sh)
            deallocate(go_sh)
         else
            deallocate(lo_nh)
            deallocate(lo_sh)
            deallocate(go_nh)
            deallocate(go_sh)
         endif
         deallocate(comb_data)

! NEIGHBOR SEARCH SHOULD HAVE SAME COMMENTS AS BILINEAR ABOVE - SEE ABOVE IF ANYTHING'S CONFUSING BELOW

      elseif(ip.eq.2) then !Neighbor search
         if(agrmet_struc(n)%gridspan.eq.1) then 
            allocate(lo_nh(agrmet_struc(n)%mo1))
            allocate(go_nh(agrmet_struc(n)%mo1))
         elseif(agrmet_struc(n)%gridspan.eq.2) then 
            allocate(lo_sh(agrmet_struc(n)%mo2))
            allocate(go_sh(agrmet_struc(n)%mo2))
         else
            allocate(lo_nh(agrmet_struc(n)%mo1))
            allocate(lo_sh(agrmet_struc(n)%mo2))
            allocate(go_nh(agrmet_struc(n)%mo1))
            allocate(go_sh(agrmet_struc(n)%mo2))
         endif
         allocate(comb_data(agrmet_struc(n)%mo1+agrmet_struc(n)%mo2))
        
         do ihemi = agrmet_struc(n)%shemi,agrmet_struc(n)%nhemi
            if(agrmet_struc(n)%global_or_hemi==1)then
               jb=int(agrmet_struc(n)%global_or_hemi*(2-ihemi)*int(jmax/2)+1)
               je=int(int(jmax)+agrmet_struc(n)%global_or_hemi*(1-ihemi)*int(jmax/2))
               jadj=agrmet_struc(n)%global_or_hemi*int(jmax/2)*(2-ihemi)
            else
               jb=1
               je=jmax
               jadj=0
            endif 
            if(agrmet_struc(n)%global_or_hemi == 0)then
               hemi=ihemi
               allocate(gi_temp(imax*jmax),li(imax*jmax))
               totdim=imax*jmax
            else
               hemi=1
               allocate(gi_temp(int(imax*jmax/2)),li(int(imax*jmax/2)))
               totdim=int(imax*jmax/2)
            endif
            li = .false.           
            gi_temp = LIS_rc%udef
            do i=1,imax
               do j=jb,je
                  if(gi(hemi,i,j).ne.udef) then 
                     jj=j-jadj
                     li(i+(jj-1)*imax) = .true.
                     gi_temp(i+(jj-1)*imax) = gi(hemi,i,j)
                  endif
               enddo
            enddo
            gridDesco(1) = LIS_rc%gridDesc(n,1) 
            gridDesco(2) = agrmet_struc(n)%hemi_nc(ihemi)
            gridDesco(3) = agrmet_struc(n)%hemi_nr(ihemi)
            gridDesco(5) = LIS_rc%gridDesc(n,5)
            gridDesco(8) = LIS_rc%gridDesc(n,8)
            gridDesco(6) = LIS_rc%gridDesc(n,6)
            gridDesco(9) = LIS_rc%gridDesc(n,9)
            gridDesco(10) = LIS_rc%gridDesc(n,10)
            if(LIS_rc%gridDesc(n,1).ne.0)gridDesco(11) = LIS_rc%gridDesc(n,11)
            !gridDesco(20) = 255
            !YDT 10/1/07
            gridDesco(20) = 0

            if(agrmet_struc(n)%gridspan.eq.1.or.agrmet_struc(n)%gridspan.eq.2) then 
               gridDesco(4) = LIS_rc%gridDesc(n,4)
               gridDesco(7) = LIS_rc%gridDesc(n,7)
            else
               if(ihemi.eq.1) then 
                  gridDesco(4) = LIS_rc%gridDesc(n,9)/2
                  gridDesco(7) = LIS_rc%gridDesc(n,7)
               else
                  gridDesco(4) = LIS_rc%gridDesc(n,4)
                  gridDesco(7) = -LIS_rc%gridDesc(n,9)/2
               endif
            endif
            if(ihemi.eq.1) then
               if(imax == agrmet_struc(n)%imax3)then 
                  call neighbor_interp(gridDesco,li,gi_temp,&
                       lo_nh,go_nh,totdim,agrmet_struc(n)%mo1,agrmet_struc(n)%rlat2_nh3, &
                       agrmet_struc(n)%rlon2_nh3,agrmet_struc(n)%n112_nh3,LIS_rc%udef,iret)
               elseif(imax == agrmet_struc(n)%imax2)then
                  call neighbor_interp(gridDesco,li,gi_temp,&
                       lo_nh,go_nh,totdim,agrmet_struc(n)%mo1,agrmet_struc(n)%rlat2_nh2, &
                       agrmet_struc(n)%rlon2_nh2,agrmet_struc(n)%n112_nh2,LIS_rc%udef,iret)
               else
                  call neighbor_interp(gridDesco,li,gi_temp,&
                       lo_nh,go_nh,totdim,agrmet_struc(n)%mo1,agrmet_struc(n)%rlat2_nh, &
                       agrmet_struc(n)%rlon2_nh,agrmet_struc(n)%n112_nh,LIS_rc%udef,iret)
               endif
            elseif(ihemi.eq.2) then 
               if(imax == agrmet_struc(n)%imax3)then
                  call neighbor_interp(gridDesco,li,gi_temp, &
                       lo_sh,go_sh,totdim,agrmet_struc(n)%mo2,agrmet_struc(n)%rlat2_sh3, & 
                       agrmet_struc(n)%rlon2_sh3,agrmet_struc(n)%n112_sh3,LIS_rc%udef,iret)
               elseif(imax == agrmet_struc(n)%imax2)then
                  call neighbor_interp(gridDesco,li,gi_temp, &
                       lo_sh,go_sh,totdim,agrmet_struc(n)%mo2,agrmet_struc(n)%rlat2_sh2, &
                       agrmet_struc(n)%rlon2_sh2,agrmet_struc(n)%n112_sh2,LIS_rc%udef,iret)
               else
                  call neighbor_interp(gridDesco,li,gi_temp, &
                       lo_sh,go_sh,totdim,agrmet_struc(n)%mo2,agrmet_struc(n)%rlat2_sh, &
                       agrmet_struc(n)%rlon2_sh,agrmet_struc(n)%n112_sh,LIS_rc%udef,iret)
               endif
            endif
            deallocate(li,gi_temp)
         enddo
         if(agrmet_struc(n)%gridspan.eq.1) then
            comb_data(1:agrmet_struc(n)%mo2) = LIS_rc%udef
            comb_data(agrmet_struc(n)%mo2+1:agrmet_struc(n)%mo1+agrmet_struc(n)%mo2) = go_nh(:)
         elseif(agrmet_struc(n)%gridspan.eq.2) then
            comb_data(1:agrmet_struc(n)%mo2) = go_sh(:)
            comb_data(agrmet_struc(n)%mo2+1:agrmet_struc(n)%mo1+agrmet_struc(n)%mo2) = LIS_rc%udef
         else
            comb_data(1:agrmet_struc(n)%mo2) = go_sh(:)
            comb_data(agrmet_struc(n)%mo2+1:agrmet_struc(n)%mo1+agrmet_struc(n)%mo2) = go_nh(:)
         endif
!just keep nhfield and shfield - low overhead, and handy for viewing in debugger
!         nhfield = RESHAPE(go_nh(1:agrmet_struc(n)%mo1),(/agrmet_struc(n)%hemi_nc(1),agrmet_struc(n)%hemi_nr(1)/))
!         shfield = RESHAPE(go_sh(1:agrmet_struc(n)%mo2),(/agrmet_struc(n)%hemi_nc(2),agrmet_struc(n)%hemi_nr(2)/))
         varfield = RESHAPE(comb_data(1:agrmet_struc(n)%mo1+agrmet_struc(n)%mo2),(/LIS_rc%lnc(n),LIS_rc%lnr(n)/))
         if(agrmet_struc(n)%gridspan.eq.1) then 
            deallocate(lo_nh)
            deallocate(go_nh)
         elseif(agrmet_struc(n)%gridspan.eq.2) then 
            deallocate(lo_sh)
            deallocate(go_sh)
         else
            deallocate(lo_nh)
            deallocate(lo_sh)
            deallocate(go_nh)
            deallocate(go_sh)
         endif
         deallocate(comb_data)
      endif
   else
      do r=LIS_rc%lnr(n)/2+1,LIS_rc%lnr(n)
         do c=1,LIS_rc%lnc(n)
            varfield(c,r) = gi(1,c,r-LIS_rc%lnr(n)/2)
         enddo
      enddo
      do r=1,LIS_rc%lnr(n)/2
         do c=1,LIS_rc%lnc(n)
            varfield(c,r) = gi(2,c,r)
         enddo
      enddo     
   endif
end subroutine interp_agrmetvar
