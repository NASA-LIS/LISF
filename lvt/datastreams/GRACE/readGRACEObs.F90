!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readGRACEObs
! \label{readGRACEObs}
!
! !INTERFACE: 
subroutine readGRACEObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,     only : LVT_rc
  use LVT_logMod
  use LVT_histDataMod
  use LVT_timeMgrMod
  use GRACE_obsMod, only : GRACEObs

#if(defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in) :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!  NOTES: 
!   The GRACE output is available at monthly intervals. So 
!   the comparisons against model data should use at least a 
!   24 hour (1day) averaging interval. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  10 Dec 2010: Sujay Kumar, Initial Specification
!  08 Jun 2018: Kristi Arsenault, Updated GRACE obs reader
! 
!EOP
  
  real        :: tws(LVT_rc%lnc, LVT_rc%lnr)
  logical     :: file_exists
  integer     :: ftn
  integer     :: timeId
  integer     :: lweId
  integer     :: tbId
  integer     :: tId
  integer     :: scaleId
  integer     :: c,r,k,t
  integer     :: mid_nc
  real        :: dt
  integer     :: currTime
  integer     :: iret
  real        :: input_data(GRACEobs(source)%nc*GRACEobs(source)%nr)
  logical*1   :: input_bitmap(GRACEobs(source)%nc*GRACEobs(source)%nr)
  real        :: output_data(LVT_rc%lnc*LVT_rc%lnr)
  logical*1   :: output_bitmap(LVT_rc%lnc*LVT_rc%lnr)
  

  tws = LVT_rc%udef

  if(GRACEobs(source)%startFlag.or.LVT_rc%resetFlag(source)) then 
     LVT_rc%resetFlag(source) = .false. 
     GRACEobs(source)%startFlag = .false. 
     
     inquire(file=trim(GRACEObs(source)%filename),exist=file_exists) 
     if(file_exists) then 
#if(defined USE_NETCDF3 || defined USE_NETCDF4) 
        write(LVT_logunit,*) '[INFO] Reading GRACE data '//&
             trim(GRACEObs(source)%filename)
        call LVT_verify(nf90_open(path=trim(GRACEObs(source)%filename),&
             mode=nf90_nowrite, ncid = ftn), &
             'nf90_open failed in readGRACEObs ')
        call LVT_verify(nf90_inq_dimid(ftn,'time',timeId),&
             'nf90_inq_dimid for time failed in readGRACEObs')
        call LVT_verify(nf90_inquire_dimension(ftn,timeId,&
             len=GRACEObs(source)%tdims),&
             'nf90_inq_dimension failed in readGRACEObs')
        
        if(allocated(GRACEObs(source)%tvals)) then 
           deallocate(GRACEObs(source)%tvals)
        endif
        if(allocated(GRACEObs(source)%time_bounds)) then 
           deallocate(GRACEObs(source)%time_bounds)
        endif
        if(allocated(GRACEObs(source)%lwe_thickness)) then 
           deallocate(GRACEObs(source)%lwe_thickness)
        endif

        allocate(GRACEObs(source)%tvals(GRACEObs(source)%tdims))
        allocate(GRACEObs(source)%time_bounds(2,GRACEObs(source)%tdims))
        allocate(GRACEObs(source)%lwe_thickness(GRACEObs(source)%nc, &
             GRACEObs(source)%nr,GRACEObs(source)%tdims))
        
        call LVT_verify(nf90_inq_varid(ftn,'lwe_thickness',lweId),&
             'nf90_inq_varid failed for lwe_thickness in readGRACEObs')
        call LVT_verify(nf90_get_var(ftn,lweId,GRACEObs(source)%lwe_thickness),&
             'nf90_get_var failed for lwe_thickness in readGRACEObs')
        
        call LVT_verify(nf90_inq_varid(ftn,'time_bounds',tbId),&
             'nf90_inq_varid failed for time_bounds in readGRACEObs')
        call LVT_verify(nf90_get_var(ftn,tbId,GRACEObs(source)%time_bounds),&
             'nf90_get_var failed for time_bounds in readGRACEObs')
        
        call LVT_verify(nf90_inq_varid(ftn,'time',tId),&
             'nf90_inq_varid failed for time in readGRACEObs')
        call LVT_verify(nf90_get_var(ftn,tId,GRACEObs(source)%tvals),&
             'nf90_get_var failed for time in readGRACEObs')
        
        call LVT_verify(nf90_close(ftn))
        
        if( trim(GRACEObs(source)%gracescalefile) .eq. "none" ) then
           write(LVT_logunit,*) '[INFO] Using unscaled GRACE data'
        else  
           write(LVT_logunit,*) '[INFO] Reading GRACE scalefactor ' &
                //trim(GRACEObs(source)%gracescalefile)
           call LVT_verify(nf90_open(path=trim(GRACEObs(source)%gracescalefile),&
                mode=nf90_nowrite, ncid = ftn), &
                'nf90_open failed in readGRACEObs scale factor file')
           if(allocated(GRACEObs(source)%scalefactor)) then 
              deallocate(GRACEObs(source)%scalefactor)
           endif

           allocate(GRACEObs(source)%scalefactor(GRACEObs(source)%nc, &
                GRACEObs(source)%nr))
           
           if( GRACEobs(source)%datasource.eq."GRACE TWS Mascon 0.5 deg") then
              call LVT_verify(nf90_inq_varid(ftn,'scale_factor',scaleId), &
                  'nf90_inq_varid failed for scale_factor in readGRACEObs')

           elseif( GRACEobs(source)%datasource.eq."GRACE TWS Original 1 deg") then
              call LVT_verify(nf90_inq_varid(ftn,'SCALE_FACTOR',scaleId), &
                  'nf90_inq_varid failed for scale_factor in readGRACEObs')
           endif
           call LVT_verify(nf90_get_var(ftn,scaleId,GRACEObs(source)%scalefactor),&
                'nf90_get_var failed for scale_factor in readGRACEObs')  
           call LVT_verify(nf90_close(ftn))
        end if
        
        write(LVT_logunit,*) '[INFO] Finished reading GRACE data '//&
             trim(GRACEObs(source)%filename)
#endif
     endif
  end if

  ! Convert GRACE 0 to 360dg domain to -180 to 180E:
  if( GRACEobs(source)%datasource.eq."GRACE TWS Mascon 0.5 deg") then
     mid_nc=360
  elseif( GRACEobs(source)%datasource.eq."GRACE TWS Original 1 deg") then
     mid_nc=180
  endif
  
  ! Estimate available time point (accounting for missing data)
  call LVT_get_julhr(LVT_rc%dyr(source),LVT_rc%dmo(source),LVT_rc%dda(source),&
       LVT_rc%dhr(source),LVT_rc%dmn(source),LVT_rc%dss(source),currTime)
  
  dt =  float((currTime-GRACEobs(source)%refTime))/24.0
  
  do k=1,GRACEobs(source)%tdims
     if(floor(GRACEobs(source)%tvals(k)).eq.dt) then 
        
        input_data = LVT_rc%udef
        
        ! No scaling applied; Apply longitude shift
        if(GRACEobs(source)%gracescalefile .eq. "none") then       
           do r=1,GRACEobs(source)%nr
              do c=1,GRACEobs(source)%nc               
                 if(GRACEobs(source)%lwe_thickness(c,r,k).ne.32767) then 
!                    if(c.le.180) then 
                    if(c.le.mid_nc) then 
!                       input_data(c+(r-1)*GRACEobs(source)%nc+180) = & 
                       input_data(c+(r-1)*GRACEobs(source)%nc+mid_nc) = & 
                            GRACEobs(source)%lwe_thickness(c,r,k)
                    else
!                       input_data(c+(r-1)*GRACEobs(source)%nc-180) = & 
                       input_data(c+(r-1)*GRACEobs(source)%nc-mid_nc) = & 
                            GRACEobs(source)%lwe_thickness(c,r,k)
                    endif
                 endif
              enddo
           enddo
            
        ! Scaling applied; Apply longitude shift
        else

           do r=1,GRACEobs(source)%nr
              do c=1,GRACEobs(source)%nc               
                 if(GRACEobs(source)%lwe_thickness(c,r,k).ne.32767) then 
!                    if(c.le.180) then 
                    if(c.le.mid_nc) then 
!                       input_data(c+(r-1)*GRACEobs(source)%nc+180) = & 
                       input_data(c+(r-1)*GRACEobs(source)%nc+mid_nc) = & 
                            (GRACEobs(source)%lwe_thickness(c,r,k) &
                            *GRACEobs(source)%scalefactor(c,r))
                    else
!                       input_data(c+(r-1)*GRACEobs(source)%nc-180) = & 
                       input_data(c+(r-1)*GRACEobs(source)%nc-mid_nc) = & 
                            (GRACEobs(source)%lwe_thickness(c,r,k) &
                            *GRACEobs(source)%scalefactor(c,r))
                    endif
                 endif
              enddo
           enddo
           
        endif
        
        input_bitmap = .false. 
        do t=1,GRACEobs(source)%nc*GRACEobs(source)%nr
           if(input_data(t).ne.LVT_rc%udef) then 
              input_bitmap(t) = .true.
           endif
        enddo
        
        ! Bilinear interpolation option:
        if( trim(GRACEObs(source)%gridtransformopt) == "bilinear" ) then
 
          call bilinear_interp(LVT_rc%gridDesc,&
             input_bitmap, input_data, &
             output_bitmap, output_data, & 
             GRACEobs(source)%nc*GRACEobs(source)%nr, & 
             LVT_rc%lnc*LVT_rc%lnr,&
             GRACEobs(source)%rlat, GRACEobs(source)%rlon, &
             GRACEobs(source)%w11, GRACEobs(source)%w12, &
             GRACEobs(source)%w21, GRACEobs(source)%w22, &
             GRACEobs(source)%n11, GRACEobs(source)%n12, &
             GRACEobs(source)%n21, GRACEobs(source)%n22, &
             LVT_rc%udef, iret)

        ! Use nearest neighbor to have product on GRACE product grid:
        elseif( trim(GRACEObs(source)%gridtransformopt) == "neighbor" ) then
          call neighbor_interp(LVT_rc%gridDesc, input_bitmap, &
             input_data, output_bitmap, output_data, &
             GRACEobs(source)%nc*GRACEobs(source)%nr, &
             LVT_rc%lnc*LVT_rc%lnr, &
             GRACEobs(source)%rlat, &
             GRACEobs(source)%rlon, &
             GRACEobs(source)%n111, &
             LVT_rc%udef, iret)
        endif
        
        ! Convert the cm to mm 
        do r=1,LVT_rc%lnr
            do c=1,LVT_rc%lnc
               if(output_data(c+(r-1)*LVT_rc%lnc).ne.LVT_rc%udef) then 
                  tws(c,r) = & 
                       output_data(c+(r-1)*LVT_rc%lnc)*10.0 !to mm. 
               else
                  tws(c,r) = LVT_rc%udef
               endif
            enddo
         enddo
         exit
      endif
   enddo

   call LVT_logSingleDataStreamVar(LVT_MOC_TWS,source,&
        tws,vlevel=1,units="mm")
   
   
 end subroutine readGRACEObs

