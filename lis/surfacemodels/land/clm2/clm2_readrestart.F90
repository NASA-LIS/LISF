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
! !ROUTINE: clm2_readrestart
! \label{clm2_readrestart}
!
! !REVISION HISTORY:
! 20 Jan 2003; Sujay Kumar Initial Specification
! 
! !INTERFACE:
subroutine clm2_readrestart()
! !USES:
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod, only : LIS_logunit, LIS_endrun
  use LIS_historyMod, only : LIS_readvar_restart
  use clm2_varpar,     only : nlevsno, nlevsoi
  use clm2_lsmMod

!
! !DESCRIPTION:
!  This program reads restart files for CLM.  This
!  includes all relevant water/energy storages and tile information. 
!
!  The routines invoked are: 
! \begin{description}
! \item[drv\_readvar\_restart](\ref{LIS_readvar_restart})
!  reads a variable from the restart file
! \end{description}
!EOP
  integer :: n
  integer :: nc,nr,ntiles
  integer :: t

  if(trim(LIS_rc%startcode).eq."restart") then 
     do n=1,LIS_rc%nnest
       
        write(LIS_logunit,*) 'Reading restart files..'
        open(40, file=clm2_struc(n)%clm_rfile,form='unformatted',status='old')
        write(LIS_logunit,*) 'CLM restart file used ',clm2_struc(n)%clm_rfile

        read(40) nc,nr,ntiles  !time, veg class, no. tiles

!------------------------------------------------------------------------
!   Check for Grid Space Conflict 
!------------------------------------------------------------------------
        if(nc.ne.LIS_rc%gnc(n).or.nr.ne.LIS_rc%gnr(n))then
           write(LIS_logunit,*)clm2_struc(n)%clm_rfile,&
                'grid space mismatch - CLM halted'
           call LIS_endrun
        endif
!------------------------------------------------------------------------
! Transfer Restart tile space to LIS tile space
!------------------------------------------------------------------------
        if(ntiles.ne.LIS_rc%glbnpatch(n,LIS_rc%lsm_index))then           
           write(LIS_logunit,*)'restart tile space mismatch, halting..'
           call LIS_endrun
        endif
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%snl)
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%frac_veg_nosno_alb)
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%h2osno)
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%h2ocan)
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%snowdp)
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%snowage)
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%frac_sno)
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%t_veg)
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%t_grnd)
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%fwet)
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%tlai)
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%tsai)
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%elai)
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%esai)
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%fsun)
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%htop)
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%hbot)
        
        do t=-nlevsno+1, 0 
           call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%dz(t))
        enddo
        do t=-nlevsno+1, 0
           call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%z(t))
        enddo
        do t=-nlevsno, 0
           call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%zi(t))
        enddo
        
        do t=-nlevsno+1, nlevsoi
           call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%t_soisno(t))
        enddo
        
        do t=1, nlevsoi
           call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%t_lake(t))
        enddo        

        do t=-nlevsno+1, nlevsoi
           call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%h2osoi_liq(t))
        enddo
        
        do t=-nlevsno+1,nlevsoi
           call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%h2osoi_ice(t))
        enddo
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%albgrd(1))
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%albgrd(2))
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%albgri(1))
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%albgri(2))
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%fabd(1))
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%fabd(2))
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%fabi(1))
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%fabi(2))
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%ftdd(1))
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%ftdd(2))
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%ftid(1))
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%ftid(2))
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%ftii(1))
        call LIS_readvar_restart(40,n,LIS_rc%lsm_index,clm2_struc(n)%clm%ftii(2))
        close(40)
     enddo
  endif

end subroutine clm2_readrestart
