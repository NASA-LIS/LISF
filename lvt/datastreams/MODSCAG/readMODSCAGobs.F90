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
! !ROUTINE: readMODSCAGobs
! \label{readMODSCAGobs}
!
! !INTERFACE:
subroutine readMODSCAGobs(source)
! 
! !USES:   
  use LVT_coreMod,    only : LVT_rc, LVT_domain
  use LVT_logMod
  use LVT_timeMgrMod, only : LVT_tick
  use LVT_histDataMod
  use map_utils
  use MODSCAG_obsMod, only : modscagobs
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
  integer,  intent(in)      :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  
!  This subroutine provides the data reader for (resampled) MODSCAG snowcover 
!  data. The routine also spatially
!  interpolates the MODSCAG data to the model (LIS) output grid 
!  and resolution. 
! 
!  NOTES: 
!   The MODSCAG output is available daily. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  09 Dec 2010: Sujay Kumar, Initial Specification
! 
!EOP

  integer, parameter           :: nc = 36000, nr=18000
  character*100                :: filename
  logical                      :: file_exists
  logical*1, allocatable           :: lb(:)
  real, allocatable                :: snfrac1(:,:)
  real                         :: snfrac(LVT_rc%lnc, LVT_rc%lnr)
  integer                      :: ftn, c,r,r1,c1
  integer                      :: ier, ivar,scfid
  real                         :: lat1, lon1
  logical                      :: alarmcheck
  real                         :: timenow

  snfrac = -9999.0
  timenow = float(LVT_rc%dhr(source))*3600 + 60*LVT_rc%dmn(source) + LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)
  if(modscagobs(source)%startflag.or.alarmCheck.or.LVT_rc%resetFlag(source)) then 
     
     LVT_rc%resetFlag(source) = .false. 
     modscagobs(source)%startflag = .false. 

     call create_modscagobs_filename(modscagobs(source)%odir,&
          LVT_rc%dyr(source), LVT_rc%dmo(source),&
          LVT_rc%dda(source), filename)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
     inquire(file=trim(filename),exist=file_exists)
     if(file_exists) then 
        write(LVT_logunit,*) '[INFO] Reading MODSCAG file ',trim(filename)
        
        allocate(lb(modscagobs(source)%modis_nc*modscagobs(source)%modis_nr))
        allocate(snfrac1(modscagobs(source)%modis_nc,modscagobs(source)%modis_nr))
        lb = .true. 

        ier = nf90_open(path=trim(filename),mode=NF90_NOWRITE, &
             ncid=ftn)
        
        if(ier.eq.0) then 
           call LVT_verify(nf90_inq_varid(ftn,'scf',scfid),&
                'nf90_inq_varid failed for scf')

           lat1 = modscagobs(source)%gridDesci(4)
           lon1 = modscagobs(source)%gridDesci(5)
           r1 = nint((lat1-0.0025)/0.005)+1
           c1 = nint((lon1-0.0025)/0.005)+1

           call LVT_verify(nf90_get_var(ftn,scfid,snfrac1, &
                start=(/c1,r1/),&
                count=(/modscagobs(source)%modis_nc, modscagobs(source)%modis_nr/)),&
                'Error in nf90_get_var: scf')
           
        endif
        call LVT_verify(nf90_close(ftn),&
             'Error in nf90_close')
#endif        
        do r=1,modscagobs(source)%modis_nr
           do c=1,modscagobs(source)%modis_nc
              if(snfrac1(c,r).lt.0.or.snfrac1(c,r).gt.100) then 
                 lb(c+(r-1)*modscagobs(source)%modis_nc) = .false.
              endif
           enddo
        enddo        

        call interp_modscag(source, &
             modscagobs(source)%modis_nc, modscagobs(source)%modis_nr, &
             snfrac1, lb, snfrac)
        

        deallocate(lb)
        deallocate(snfrac1)
     else
        snfrac = -9999
     endif
  else
     snfrac = -9999.0
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_snowcover,source,snfrac,vlevel=1,units="-")

end subroutine readMODSCAGobs

!BOP
! 
! !ROUTINE: interp_modscag
!  \label{interp_modscag}
!
! !INTERFACE: 
subroutine interp_modscag(source, nc, nr, var_input, lb, var_output)
! 
! !USES:   
  use LVT_coreMod
  use MODSCAG_obsMod, only : modscagobs

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!   This subroutine spatially interpolates the MODSCAG variable to the
!   model (LIS) output grid and resolution, using a bilinear interpolation
!   approach. 
! 
!   The arguments are: 
!   \begin{description}
!    \item[nc]      number of columns in the input (MODSCAG) grid
!    \item[nr]      number of rows in the input (MODSCAG) grid
!    \item[var_input] input variable to be interpolated
!    \item[lb]        input bitmap (true//false)
!    \item[var_output] resulting interplated field
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS:   
  integer            :: source
  integer            :: nc
  integer            :: nr
  real               :: var_input(nc*nr)
  logical*1          :: lb(nc*nr)
  real               :: var_output(LVT_rc%lnc, LVT_rc%lnr)
!EOP
  integer            :: iret
  integer            :: c,r
  logical*1          :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real               :: go(LVT_rc%lnc*LVT_rc%lnr)

  var_output = LVT_rc%udef
  if(LVT_isAtAfinerResolution(0.005)) then
     call bilinear_interp(LVT_rc%gridDesc,lb, var_input,&
          lo,go,nc*nr,LVT_rc%lnc*LVT_rc%lnr,&
          modscagobs(source)%rlat,modscagobs(source)%rlon,&
          modscagobs(source)%w11,modscagobs(source)%w12,&
          modscagobs(source)%w21,modscagobs(source)%w22,&
          modscagobs(source)%n11,modscagobs(source)%n12,&
          modscagobs(source)%n21,modscagobs(source)%n22,LVT_rc%udef,iret)
     
  else
     call upscaleByAveraging(&
          modscagobs(source)%modis_nc*modscagobs(source)%modis_nr, &
          LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
          modscagobs(source)%n11, lb, &
          var_input, lo, go)
  endif

  do r=1, LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(go(c+(r-1)*LVT_rc%lnc).ne.-9999.0) then 
           var_output(c,r) = go(c+(r-1)*LVT_rc%lnc)/100.0
        endif
     enddo
  enddo
  
end subroutine interp_modscag

!BOP
! 
! !ROUTINE: create_modscagobs_filename
! \label{create_modscagobs_filename}
!
! !INTERFACE: 
subroutine create_modscagobs_filename(odir, yr, mo, da, filename)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine creates the MODSCAG filename based on the given 
!  date (year, month, day and hour)
! 
!  The arguments are: 
!  \begin{description}
!   \item[odir]      MODSCAG base directory
!   \item[yr]        year of data
!   \item[mo]        month of data
!   \item[da]        day of data
!   \item[filename]  Name of the MODSCAG file
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS:  
  character(len=*)  :: odir
  integer           :: yr
  integer           :: mo
  integer           :: da
  character(len=*)  :: filename
! 
!EOP

  character*4       :: fyr
  character*2       :: fmo
  character*2       :: fda

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da

  filename = trim(odir)//'/'//trim(fyr)//'/griddata/MOD09GA.'//&
       trim(fyr)//trim(fmo)//trim(fda)//&
       '.lat_0-90_lon_0-180.snow_fraction_5e-3deg.nc'

end subroutine create_modscagobs_filename
  
