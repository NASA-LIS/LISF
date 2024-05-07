!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: jules5x_setup
! \label{jules5x_setup}
!
! !REVISION HISTORY:
! 16 May 2016; Shugong Wang; initial implementation for JULES 4.3
! 01 Feb 2018; Shugong Wang; updated for JULES 5.0 
! 28 Nov 2018; Shugong Wang; updated for JULES 5.2 
! 12 Dec 2018; Shugong Wang; updated for JULES 5.3 
! 08 Jul 2019; Shugong Wang; updated for JULES 5.6
! 
! !INTERFACE:
subroutine jules5x_setup()
! !USES:
  use LIS_logMod,    only: LIS_logunit, LIS_verify, LIS_endrun
  use LIS_fileIOMod, only: LIS_read_param
  use LIS_coreMod,   only: LIS_rc, LIS_surface
  use jules_surface_mod,      only: l_aggregate
  use jules5x_lsmMod
! !ARGUMENTS: 
!
! !DESCRIPTION:
!
!  Complete the setup routines for Template option (forcing-only) 
! 
!EOP

  implicit none
  integer :: n
  integer           :: mtype
  integer           :: t, k
  integer           :: col, row
  real, allocatable :: placeholder(:,:)
  
  mtype = LIS_rc%lsm_index
  
  do n=1, LIS_rc%nnest
      ! allocate memory for place holder for #n nest
      allocate(placeholder(LIS_rc%lnc(n), LIS_rc%lnr(n)))
      
      !------------------------------------!
      ! reading spatial spatial parameters !
      !------------------------------------!
      ! vegetype takes value from the LIS built-in parameter vegt
      write(LIS_logunit,*) "JULES: retrieve parameter VEGETYPE from LIS"
      do t=1, LIS_rc%npatch(n, mtype)
        if(l_aggregate) then
          jules5x_struc(n)%jules5x(t)%pft = 1
        else
          jules5x_struc(n)%jules5x(t)%pft = LIS_surface(n, mtype)%tile(t)%vegt
        endif
      enddo

      write(LIS_logunit,*) "JULES: reading parameter B from ", trim(LIS_rc%paramfile(n))
      call LIS_read_param(n, "JULES_B", placeholder)
      do t = 1, LIS_rc%npatch(n, mtype)
          col = LIS_surface(n, mtype)%tile(t)%col
          row = LIS_surface(n, mtype)%tile(t)%row
          jules5x_struc(n)%jules5x(t)%b = placeholder(col, row)
      enddo 
      
      write(LIS_logunit,*) "JULES: reading parameter SATHH from ", trim(LIS_rc%paramfile(n))
      call LIS_read_param(n, "JULES_SATHH", placeholder)
      do t = 1, LIS_rc%npatch(n, mtype)
          col = LIS_surface(n, mtype)%tile(t)%col
          row = LIS_surface(n, mtype)%tile(t)%row
          jules5x_struc(n)%jules5x(t)%sathh = placeholder(col, row)
          jules5x_struc(n)%jules5x(t)%p_s_sathh(:) = jules5x_struc(n)%jules5x(t)%sathh 
      enddo 
      
      write(LIS_logunit,*) "JULES: reading parameter SATCON from ", trim(LIS_rc%paramfile(n))
      call LIS_read_param(n, "JULES_SATCON", placeholder)
      do t = 1, LIS_rc%npatch(n, mtype)
          col = LIS_surface(n, mtype)%tile(t)%col
          row = LIS_surface(n, mtype)%tile(t)%row
          jules5x_struc(n)%jules5x(t)%satcon = placeholder(col, row)
      enddo 
      

      write(LIS_logunit,*) "JULES: reading parameter SM_SAT from ", trim(LIS_rc%paramfile(n))
      call LIS_read_param(n, "JULES_SM_SAT", placeholder)
      do t = 1, LIS_rc%npatch(n, mtype)
          col = LIS_surface(n, mtype)%tile(t)%col
          row = LIS_surface(n, mtype)%tile(t)%row
          jules5x_struc(n)%jules5x(t)%sm_sat = placeholder(col, row)
      enddo 

      write(LIS_logunit,*) "JULES: reading parameter SM_CRIT from ", trim(LIS_rc%paramfile(n))
      call LIS_read_param(n, "JULES_SM_CRIT", placeholder)
      do t = 1, LIS_rc%npatch(n, mtype)
          col = LIS_surface(n, mtype)%tile(t)%col
          row = LIS_surface(n, mtype)%tile(t)%row
          jules5x_struc(n)%jules5x(t)%sm_crit = placeholder(col, row)
      enddo 

      write(LIS_logunit,*) "JULES: reading parameter SM_WILT from ", trim(LIS_rc%paramfile(n))
      call LIS_read_param(n, "JULES_SM_WILT", placeholder)
      do t = 1, LIS_rc%npatch(n, mtype)
          col = LIS_surface(n, mtype)%tile(t)%col
          row = LIS_surface(n, mtype)%tile(t)%row
          jules5x_struc(n)%jules5x(t)%sm_wilt = placeholder(col, row)
      enddo 

      write(LIS_logunit,*) "JULES: reading parameter HCAP from ", trim(LIS_rc%paramfile(n))
      call LIS_read_param(n, "JULES_HCAP", placeholder)
      do t = 1, LIS_rc%npatch(n, mtype)
          col = LIS_surface(n, mtype)%tile(t)%col
          row = LIS_surface(n, mtype)%tile(t)%row
          jules5x_struc(n)%jules5x(t)%hcap = placeholder(col, row)
      enddo 

      write(LIS_logunit,*) "JULES: reading parameter HCON from ", trim(LIS_rc%paramfile(n))
      call LIS_read_param(n, "JULES_HCON", placeholder)
      do t = 1, LIS_rc%npatch(n, mtype)
          col = LIS_surface(n, mtype)%tile(t)%col
          row = LIS_surface(n, mtype)%tile(t)%row
          jules5x_struc(n)%jules5x(t)%hcon = placeholder(col, row)
      enddo 

      write(LIS_logunit,*) "JULES: reading parameter ALBSOIL from ", trim(LIS_rc%paramfile(n))
      call LIS_read_param(n, "JULES_ALBSOIL", placeholder)
      do t = 1, LIS_rc%npatch(n, mtype)
          col = LIS_surface(n, mtype)%tile(t)%col
          row = LIS_surface(n, mtype)%tile(t)%row
          jules5x_struc(n)%jules5x(t)%albsoil = placeholder(col, row)
      enddo 
      
      write(LIS_logunit,*) "JULES: reading parameter TOPMODEL ti_mean from ", trim(LIS_rc%paramfile(n))
      call LIS_read_param(n, "JULES_TI_MEAN", placeholder)
      do t = 1, LIS_rc%npatch(n, mtype)
          col = LIS_surface(n, mtype)%tile(t)%col
          row = LIS_surface(n, mtype)%tile(t)%row
          jules5x_struc(n)%jules5x(t)%ti_mean = placeholder(col, row)
      enddo 
      
      write(LIS_logunit,*) "JULES: reading parameter TOPMODEL ti_sig from ", trim(LIS_rc%paramfile(n))
      call LIS_read_param(n, "JULES_TI_SIG", placeholder)
      do t = 1, LIS_rc%npatch(n, mtype)
          col = LIS_surface(n, mtype)%tile(t)%col
          row = LIS_surface(n, mtype)%tile(t)%row
          jules5x_struc(n)%jules5x(t)%ti_sig = placeholder(col, row)
          if (jules5x_struc(n)%jules5x(t)%ti_sig ==0) then
            jules5x_struc(n)%jules5x(t)%ti_sig = 0.1
          endif
      enddo 
      
      write(LIS_logunit,*) "JULES: reading parameter TOPMODEL fexp from ", trim(LIS_rc%paramfile(n))
      call LIS_read_param(n, "JULES_FEXP", placeholder)
      do t = 1, LIS_rc%npatch(n, mtype)
          col = LIS_surface(n, mtype)%tile(t)%col
          row = LIS_surface(n, mtype)%tile(t)%row
          jules5x_struc(n)%jules5x(t)%fexp = placeholder(col, row)
      enddo 
        
      if(l_aggregate) then 
        write(LIS_logunit,*) "JULES: reading LANDCOVER from ", trim(LIS_rc%paramfile(n))
        do k = 1, jules5x_struc(n)%ntype
            call jules5x_read_multilevel_param(n, "LANDCOVER", k, placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                jules5x_struc(n)%jules5x(t)%surft_frac(k) = placeholder(col, row)
            enddo 
        enddo 
      endif 
  enddo

end subroutine jules5x_setup
