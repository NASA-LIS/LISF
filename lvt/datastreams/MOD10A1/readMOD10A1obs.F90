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
! !ROUTINE: readMOD10A1obs
! \label{readMOD10A1obs}
!
! !INTERFACE:
subroutine readMOD10A1obs(source)
! 
! !USES:   
  use LVT_coreMod,    only : LVT_rc, LVT_domain
  use LVT_logMod
  use LVT_timeMgrMod, only : LVT_tick
  use LVT_histDataMod
  use map_utils
  use MOD10A1_obsMod, only : mod10a1obs
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
!  This subroutine provides the data reader for (resampled) MOD10A1 snowcover 
!  data. The routine also spatially
!  interpolates the MOD10A1 data to the model (LIS) output grid 
!  and resolution. 
! 
!  NOTES: 
!   The MOD10A1 output is available daily. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  09 Dec 2010: Sujay Kumar, Initial Specification
! 
!EOP

  integer, parameter           :: nc = 36000, nr=15000
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
  if(mod10a1obs(source)%startflag.or.alarmCheck.or.LVT_rc%resetFlag(source)) then 
     
     LVT_rc%resetFlag(source) = .false. 
     mod10a1obs(source)%startflag = .false. 

     call create_mod10a1obs_filename(mod10a1obs(source)%odir,&
          LVT_rc%dyr(source), LVT_rc%dmo(source),&
          LVT_rc%dda(source), filename)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
     inquire(file=trim(filename),exist=file_exists)
     if(file_exists) then 
        write(LVT_logunit,*) '[INFO] Reading MOD10A1 file ',trim(filename)
        
        allocate(lb(mod10a1obs(source)%modis_nc*mod10a1obs(source)%modis_nr))
        allocate(snfrac1(mod10a1obs(source)%modis_nc,mod10a1obs(source)%modis_nr))
        lb = .true. 

        ier = nf90_open(path=trim(filename),mode=NF90_NOWRITE, &
             ncid=ftn)
        
        if(ier.eq.0) then 
           call LVT_verify(nf90_inq_varid(ftn,'scf',scfid),&
                'nf90_inq_varid failed for scf')

           lat1 = mod10a1obs(source)%gridDesci(4)
           lon1 = mod10a1obs(source)%gridDesci(5)
           r1 = nint((lat1+59.995)/0.01)+1
           c1 = nint((lon1+179.995)/0.01)+1

           call LVT_verify(nf90_get_var(ftn,scfid,snfrac1, &
                   start=(/c1,r1/),&
                   count=(/mod10a1obs(source)%modis_nc, mod10a1obs(source)%modis_nr/)),&
                   'Error in nf90_get_var: scf')
              
        endif
        call LVT_verify(nf90_close(ftn),&
             'Error in nf90_close')
#endif        
        do r=1,mod10a1obs(source)%modis_nr
           do c=1,mod10a1obs(source)%modis_nc
              if(snfrac1(c,r).lt.0.or.snfrac1(c,r).gt.100) then 
                 lb(c+(r-1)*mod10a1obs(source)%modis_nc) = .false.
              endif
           enddo
        enddo        

        call interp_mod10a1(source, &
             mod10a1obs(source)%modis_nc, mod10a1obs(source)%modis_nr, &
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

end subroutine readMOD10A1obs

!BOP
! 
! !ROUTINE: interp_mod10a1
!  \label{interp_mod10a1}
!
! !INTERFACE: 
subroutine interp_mod10a1(source, nc, nr, var_input, lb, var_output)
! 
! !USES:   
  use LVT_coreMod
  use MOD10A1_obsMod, only : mod10a1obs

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!   This subroutine spatially interpolates the MOD10A1 variable to the
!   model (LIS) output grid and resolution, using a bilinear interpolation
!   approach. 
! 
!   The arguments are: 
!   \begin{description}
!    \item[nc]      number of columns in the input (MOD10A1) grid
!    \item[nr]      number of rows in the input (MOD10A1) grid
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
  if(LVT_isAtAfinerResolution(0.01)) then
     call bilinear_interp(LVT_rc%gridDesc,lb, var_input,&
          lo,go,nc*nr,LVT_rc%lnc*LVT_rc%lnr,&
          mod10a1obs(source)%rlat,mod10a1obs(source)%rlon,&
          mod10a1obs(source)%w11,mod10a1obs(source)%w12,&
          mod10a1obs(source)%w21,mod10a1obs(source)%w22,&
          mod10a1obs(source)%n11,mod10a1obs(source)%n12,&
          mod10a1obs(source)%n21,mod10a1obs(source)%n22,LVT_rc%udef,iret)
     
  else
     call upscaleByAveraging(&
          mod10a1obs(source)%modis_nc*mod10a1obs(source)%modis_nr, &
          LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
          mod10a1obs(source)%n11, lb, &
          var_input, lo, go)
  endif

  do r=1, LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(go(c+(r-1)*LVT_rc%lnc).ne.-9999.0) then 
           var_output(c,r) = go(c+(r-1)*LVT_rc%lnc)/100.0
        endif
     enddo
  enddo
  
end subroutine interp_mod10a1

!BOP
! 
! !ROUTINE: create_mod10a1obs_filename
! \label{create_mod10a1obs_filename}
!
! !INTERFACE: 
subroutine create_mod10a1obs_filename(odir, yr, mo, da, filename)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine creates the MOD10A1 filename based on the given 
!  date (year, month, day and hour)
! 
!  The arguments are: 
!  \begin{description}
!   \item[odir]      MOD10A1 base directory
!   \item[yr]        year of data
!   \item[mo]        month of data
!   \item[da]        day of data
!   \item[filename]  Name of the MOD10A1 file
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

  filename = trim(odir)//'/'//trim(fyr)//'/MOD10A1_'//trim(fyr)//trim(fmo)//trim(fda)//&
       '_c5_1km.nc4'

end subroutine create_mod10a1obs_filename
  
