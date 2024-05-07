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
! !ROUTINE: interp_agrradpsvar 
!  \label{interp_agrradpsvar} 
! 
! !INTERFACE:
subroutine interp_agrradpsvar(n,ip,gi,udef,varfield)
! !USES:
  use LIS_coreMod,         only : LIS_rc
  use agrradps_forcingMod, only : agrradps_struc
  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n 
  integer, intent(in)    :: ip
  real, intent(in)       :: gi(2,agrradps_struc(n)%imax,agrradps_struc(n)%jmax)
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
!  \begin{description}
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[neighbor\_interp](\ref{neighbor_interp}) \newline
!    spatially interpolate the forcing data using neighbor interpolation
! \end{description}
!EOP
  real                   :: gi_temp(agrradps_struc(n)%mi)
  logical*1              :: li(agrradps_struc(n)%mi)
  integer                :: ihemi,i,j,iret
  real                   :: gridDesco(200)
  logical*1, allocatable :: lo_nh(:),lo_sh(:)
  real, allocatable      :: go_nh(:),go_sh(:)
  real, allocatable      :: comb_data(:)
  integer                :: c,r

  if(agrradps_struc(n)%interp.gt.0) then 
     if(ip.eq.1) then 

        if(agrradps_struc(n)%gridspan.eq.1) then 
           allocate(lo_nh(agrradps_struc(n)%mo1))
           allocate(go_nh(agrradps_struc(n)%mo1))
        elseif(agrradps_struc(n)%gridspan.eq.2) then 
           allocate(lo_sh(agrradps_struc(n)%mo2))
           allocate(go_sh(agrradps_struc(n)%mo2))
        else
           allocate(lo_nh(agrradps_struc(n)%mo1))
           allocate(lo_sh(agrradps_struc(n)%mo2))
           allocate(go_nh(agrradps_struc(n)%mo1))
           allocate(go_sh(agrradps_struc(n)%mo2))
        endif
      allocate(comb_data(1:agrradps_struc(n)%mo1+agrradps_struc(n)%mo2))
        
        do ihemi = agrradps_struc(n)%shemi,agrradps_struc(n)%nhemi
           li = .false.
           do i=1,agrradps_struc(n)%imax
              do j=1,agrradps_struc(n)%jmax
                 if(gi(ihemi,i,j).ne.udef) then 
                    li(i+(j-1)*agrradps_struc(n)%imax) = .true.
                    gi_temp(i+(j-1)*agrradps_struc(n)%imax) = gi(ihemi,i,j)
                 endif
              enddo
           enddo
           gridDesco(1) = 0 
           gridDesco(2) = agrradps_struc(n)%hemi_nc(ihemi)
           gridDesco(3) = agrradps_struc(n)%hemi_nr(ihemi)
           gridDesco(5) = LIS_rc%gridDesc(n,5)
           gridDesco(8) = LIS_rc%gridDesc(n,8)
           gridDesco(6) = LIS_rc%gridDesc(n,6)
           gridDesco(9) = LIS_rc%gridDesc(n,9)
           gridDesco(10) = LIS_rc%gridDesc(n,10)
           !YDT 10/2/07
           !gridDesco(20) = 255
           gridDesco(20) = 0 
           if(agrradps_struc(n)%gridspan.eq.1 .or. &
              agrradps_struc(n)%gridspan.eq.2) then 
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
              call bilinear_interp(gridDesco,li,gi_temp,&
                   lo_nh,go_nh,agrradps_struc(n)%mi,agrradps_struc(n)%mo1,&
                   agrradps_struc(n)%rlat1_nh,agrradps_struc(n)%rlon1_nh,&
                   agrradps_struc(n)%w111_nh,agrradps_struc(n)%w121_nh,&
                   agrradps_struc(n)%w211_nh,agrradps_struc(n)%w221_nh,&
                   agrradps_struc(n)%n111_nh,agrradps_struc(n)%n121_nh,&
                   agrradps_struc(n)%n211_nh,agrradps_struc(n)%n221_nh,LIS_rc%udef,&
                   iret)
           elseif(ihemi.eq.2) then 
              call bilinear_interp(gridDesco,li,gi_temp,&
                   lo_sh,go_sh,agrradps_struc(n)%mi,agrradps_struc(n)%mo2,&
                   agrradps_struc(n)%rlat1_sh, agrradps_struc(n)%rlon1_sh,&
                   agrradps_struc(n)%w111_sh,agrradps_struc(n)%w121_sh,&
                   agrradps_struc(n)%w211_sh,agrradps_struc(n)%w221_sh,&
                   agrradps_struc(n)%n111_sh,agrradps_struc(n)%n121_sh,&
                   agrradps_struc(n)%n211_sh,agrradps_struc(n)%n221_sh,&
                   LIS_rc%udef,iret)
           endif
        end do
        if(agrradps_struc(n)%gridspan.eq.1) then
           comb_data(1:agrradps_struc(n)%mo1) = go_nh(:)
        elseif(agrradps_struc(n)%gridspan.eq.2) then 
           comb_data(1:agrradps_struc(n)%mo2) = go_sh(:)
        else
           comb_data(1:agrradps_struc(n)%mo2) = go_sh(:)
           comb_data(agrradps_struc(n)%mo2+1:agrradps_struc(n)%mo1+agrradps_struc(n)%mo2) = go_nh(:)
        endif

        do r=1,LIS_rc%lnr(n)
         do c=1,LIS_rc%lnc(n)
            varfield(c,r) = comb_data(c+(r-1)*LIS_rc%lnc(n))
         enddo
        enddo
        
        if(agrradps_struc(n)%gridspan.eq.1) then 
           deallocate(lo_nh)
           deallocate(go_nh)
        elseif(agrradps_struc(n)%gridspan.eq.2) then 
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

end subroutine interp_agrradpsvar

