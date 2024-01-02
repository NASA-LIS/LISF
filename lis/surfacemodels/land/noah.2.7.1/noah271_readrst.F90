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
! !ROUTINE: noah271_readrst
! \label{noah271_readrst}
!
! !REVISION HISTORY:
!  05 Sep 2001: Brian Cosgrove; Modified code to use Dag Lohmann's NOAA
!               initial conditions if necessary.  This is controlled with
!               local variable NOAAIC.  Normally set to 0 in this subroutine
!               but set to 1 if want to use Dag's NOAA IC's.  Changed output
!               directory structure, and commented out if-then check so that
!               directory is always made.
!  28 Apr 2002: Kristi Arsenault; Added Noah2.5 LSM into LDAS
!  28 May 2002: Kristi Arsenault; For STARTCODE=4, corrected SNEQV values  
!                and put SMC, SH2O, STC limit for GDAS and GEOS forcing.
!  27 Oct 2010: David Mocko, changes for Noah2.7.1 in LIS6.1
!
! !INTERFACE:
subroutine noah271_readrst
! !USES:
  use LIS_coreMod
  use LIS_historyMod, only : LIS_readvar_restart
  use LIS_logMod,     only : LIS_logunit, LIS_endrun, &
       LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify
  use noah271_lsmMod
  use LIS_tbotAdjustMod, only: LIS_readTmnUpdateRestart
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
!
! !DESCRIPTION:
!  This program reads restart files for Noah2.7.1.  This
!  includes all relevant water/energy storages and tile information. 
!  The following is the list of variables specified in the Noah2.7.1 
!  restart file: 
!  \begin{verbatim}
!   nc,nr,ntiles    - grid and tile space dimensions 
!   t1           - Noah2.7.1 skin temperature
!   cmc          - Noah2.7.1 canopy moisture storage
!   snowh        - Noah2.7.1 snow depth
!   sneqv        - Noah2.7.1 snow water equivalent
!   stc          - Noah2.7.1 soil temperature (for each layer)
!   smc          - Noah2.7.1 soil moisture (for each layer) 
!   sh2o         - Noah2.7.1 liquid only soil moisture (for each layer)
!   ch           - Noah2.7.1 heat/moisture exchange coefficient
!   cm           - Noah2.7.1 momentum exchange coefficient
!  \end{verbatim}
!
!  The routines invoked are: 
! \begin{description}
! \item[LIS\_readvar\_restart](\ref{LIS_readvar_restart}) \newline
!  reads a variable from the restart file
! \item[noah\_coldstart](\ref{noah271_coldstart}) \newline
!   initializes the Noah2.7.1 state variables
! \item[noah271\_read\_paramrst](\ref{noah271_read_paramrst}) \newline
!   reads the OPTUE output to read the optimized LSM parameters 
! \end{description}
!EOP
  implicit none      
  integer           :: t,l
  integer           :: nc,nr,ntiles
  integer           :: n
  integer           :: ftn
  integer           :: status
  real, allocatable :: tmptilen(:)
  logical           :: file_exists
  character*20      :: wformat

  wformat = "netcdf"
  if(trim(LIS_rc%startcode).eq."restart") then 
     do n=1,LIS_rc%nnest
        allocate(tmptilen(LIS_rc%npatch(n,LIS_rc%lsm_index)))
!-------------------------------------------------------------------------
! Read Active Archive File
!-------------------------------------------------------------------------
        inquire(file=trim(noah271_struc(n)%rfile),exist=file_exists) 

        if(.not.file_exists) then 
           write(LIS_logunit,*) 'Noah2.7.1 restart file ',             &
                         trim(noah271_struc(n)%rfile),' does not exist '
           write(LIS_logunit,*) 'Program stopping ...'
           call LIS_endrun()
           
        endif
        write(LIS_logunit,*)                        &
                  'Noah2.7.1 restart file used: ',noah271_struc(n)%rfile
        if(trim(wformat).eq."binary") then 
           ftn = LIS_getNextUnitNumber()
           open(ftn,file=noah271_struc(n)%rfile,form='unformatted')        
           read(ftn) nc,nr,ntiles  !time, veg class, no. tiles

!------------------------------------------------------------------------
!   Check for Grid Space Conflict 
!------------------------------------------------------------------------
           if(nc.ne.LIS_rc%gnc(n) .or. nr.ne.LIS_rc%gnr(n))then
              write(LIS_logunit,*) noah271_struc(n)%rfile,                &
                   'grid space mismatch - Noah2.7.1 halted'
              call LIS_endrun
           endif
!------------------------------------------------------------------------
! Transfer Restart tile space to LIS tile space
!------------------------------------------------------------------------
           if(ntiles.ne.LIS_rc%glbntiles_red(n))then           
              write(LIS_logunit,*) 'restart tile space mismatch, halting..'
              call LIS_endrun
           endif
        elseif(trim(wformat).eq."netcdf") then 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
           status=nf90_open(path=trim(noah271_struc(n)%rfile),&
                mode=NF90_NOWRITE,ncid=ftn)
           call LIS_verify(status,&
                'Error opening file '//trim(noah271_struc(n)%rfile))
#endif
        endif

        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,&
             noah271_struc(n)%noah%t1, &
             varname="T1",wformat=trim(wformat))
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,&
             noah271_struc(n)%noah%cmc,&
             varname="CMC", wformat=trim(wformat))
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,&
             noah271_struc(n)%noah%snowh,&
             varname="SNOWH",wformat=trim(wformat))
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,&
             noah271_struc(n)%noah%sneqv,&
             varname="SNEQV",wformat=trim(wformat))

        do l=1,noah271_struc(n)%nslay
           call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,&
                tmptilen, varname="STC",&
                dim=l,vlevels = noah271_struc(n)%nslay, &
                wformat=trim(wformat))
           do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
              noah271_struc(n)%noah(t)%stc(l) = tmptilen(t)
           enddo
        enddo

        do l=1,noah271_struc(n)%nslay
           call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,&
                tmptilen, varname="SMC",&
                dim=l,vlevels = noah271_struc(n)%nslay, &
                wformat=trim(wformat))
           do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
              noah271_struc(n)%noah(t)%smc(l) = tmptilen(t)
           enddo
        enddo

        do l=1,noah271_struc(n)%nslay
           call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,&
                tmptilen, varname="SH2O",&
                dim=l,vlevels = noah271_struc(n)%nslay, &
                wformat=trim(wformat))
           do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
              if(tmptilen(t).lt.0.0001)  tmptilen(t) = 0.0001
              noah271_struc(n)%noah(t)%sh2o(l) = tmptilen(t)
           enddo
        enddo
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,&
             noah271_struc(n)%noah%ch,&
             varname="CH",wformat=trim(wformat))
        call LIS_readvar_restart(ftn,n, LIS_rc%lsm_index,&
             noah271_struc(n)%noah%cm,&
             varname="CM",wformat=trim(wformat))

        !Read in variables for dynamic deep soil temperature
        call LIS_readTmnUpdateRestart(n,ftn,wformat)

        if(trim(wformat).eq."binary") then 
           call LIS_releaseUnitNumber(ftn)        
        elseif(trim(wformat).eq."netcdf") then 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
           status = nf90_close(ftn)
           call LIS_verify(status,'Error in nf90_close in noah271_readrst')
#endif     
        endif
        deallocate(tmptilen)

     enddo
  elseif(trim(LIS_rc%startcode).eq."coldstart") then 
     call noah271_coldstart()
  endif

  do n=1,LIS_rc%nnest
     noah271_struc(n)%count = 0 
  enddo

  do n=1,LIS_rc%nnest
     if(noah271_struc(n)%param_rst.eq.1) then 
        call noah271_read_paramrst(n)        
     endif
  enddo

end subroutine noah271_readrst

!BOP
! !ROUTINE: noah271_read_paramrst
! \label{noah271_read_paramrst}
!
! !INTERFACE: 
subroutine noah271_read_paramrst(n)
! !USES:   
  use LIS_coreMod
  use LIS_logMod,       only : LIS_getNextUnitNumber, LIS_releaseUnitNumber
  use LIS_historyMod,   only : LIS_readvar_gridded
  use noah271_lsmMod

  implicit none
! !ARGUMENTS: 
  integer,             intent(in)     :: n 
!
! !DESCRIPTION: 
!  This subroutine reads the output file of the OPT/UE algorithm to read
!  the optimized set of model parameters. The optimized parameters are
!  then used to overwrite the LSM's default parameters. 
!
!EOP
  integer                             :: ftn
  integer                             :: nparam
  real,  allocatable                      :: params(:,:,:)
  integer                             :: t,k,col,row
  character*100                       :: param_name

  ftn = LIS_getNextUnitNumber()
  open(ftn,file=trim(noah271_struc(n)%prstfile),form='unformatted') 
  
  read(ftn) nparam
  allocate (params(nparam, LIS_rc%lnc(n),LIS_rc%lnr(n)))
  
  do k=1,nparam
     read(ftn) param_name
     call LIS_readvar_gridded(ftn,n,params(k,:,:))
     if(trim(param_name).eq."SMCMAX") then 
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
           row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
           noah271_struc(n)%noah(t)%smcmax = params(k,col,row)
        enddo
     endif
     if(trim(param_name).eq."PSISAT") then 
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
           row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
           noah271_struc(n)%noah(t)%psisat = params(k,col,row)
        enddo
     endif
     if(trim(param_name).eq."DKSAT") then 
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
           row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
           noah271_struc(n)%noah(t)%dksat = params(k,col,row)
        enddo
     endif
     if(trim(param_name).eq."DWSAT") then 
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
           row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
           noah271_struc(n)%noah(t)%dwsat = params(k,col,row)
        enddo
     endif
     if(trim(param_name).eq."BEXP") then 
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
           row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
           noah271_struc(n)%noah(t)%bexp = params(k,col,row)
        enddo
     endif
     if(trim(param_name).eq."QUARTZ") then 
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
           row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
           noah271_struc(n)%noah(t)%quartz = params(k,col,row)
        enddo
     endif
     if(trim(param_name).eq."RSMIN") then 
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
           row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
           noah271_struc(n)%noah(t)%rsmin = params(k,col,row)
        enddo
     endif
     if(trim(param_name).eq."RGL") then 
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
           row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
           noah271_struc(n)%noah(t)%rgl = params(k,col,row)
        enddo
     endif
     if(trim(param_name).eq."HS") then 
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
           row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
           noah271_struc(n)%noah(t)%hs = params(k,col,row)
        enddo
     endif
     if(trim(param_name).eq."Z0") then 
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
           row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
           noah271_struc(n)%noah(t)%z0 = params(k,col,row)
        enddo
     endif
     if(trim(param_name).eq."LAI") then 
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
           row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
           noah271_struc(n)%noah(t)%lai = params(k,col,row)
        enddo
     endif
     if(trim(param_name).eq."CFACTR") then 
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
           row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
           noah271_struc(n)%noah(t)%cfactr = params(k,col,row)
        enddo
     endif
     if(trim(param_name).eq."CMCMAX") then 
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
           row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
           noah271_struc(n)%noah(t)%cmcmax = params(k,col,row)
        enddo
     endif
     if(trim(param_name).eq."SBETA") then 
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
           row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
           noah271_struc(n)%noah(t)%sbeta = params(k,col,row)
        enddo
     endif
     if(trim(param_name).eq."RSMAX") then 
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
           row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
           noah271_struc(n)%noah(t)%rsmax = params(k,col,row)
        enddo
     endif
     if(trim(param_name).eq."TOPT") then 
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
           row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
           noah271_struc(n)%noah(t)%topt = params(k,col,row)
        enddo
     endif
     if(trim(param_name).eq."REFDK") then 
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
           row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
           noah271_struc(n)%noah(t)%refdk = params(k,col,row)
        enddo
     endif
     if(trim(param_name).eq."FXEXP") then 
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
           row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
           noah271_struc(n)%noah(t)%fxexp = params(k,col,row)
        enddo
     endif
     if(trim(param_name).eq."REFKDT") then 
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
           row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
           noah271_struc(n)%noah(t)%refkdt = params(k,col,row)
        enddo
     endif
     if(trim(param_name).eq."CZIL") then 
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
           row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
           noah271_struc(n)%noah(t)%czil = params(k,col,row)
        enddo
     endif
     if(trim(param_name).eq."CSOIL") then 
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
           row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
           noah271_struc(n)%noah(t)%csoil = params(k,col,row)
        enddo
     endif
     if(trim(param_name).eq."FRZK") then 
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
           row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
           noah271_struc(n)%noah(t)%frzk = params(k,col,row)
        enddo
     endif
     if(trim(param_name).eq."SNUP") then 
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
           row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
           noah271_struc(n)%noah(t)%snup = params(k,col,row)
        enddo
     endif
     if(trim(param_name).eq."SMCREF") then 
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
           row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
           noah271_struc(n)%noah(t)%smcref = params(k,col,row)
        enddo
     endif
     if(trim(param_name).eq."SMCWLT") then 
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
           row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
           noah271_struc(n)%noah(t)%smcwlt = params(k,col,row)
        enddo
     endif

  enddo

  deallocate(params)
  call LIS_releaseUnitNumber(ftn)

end subroutine noah271_read_paramrst
