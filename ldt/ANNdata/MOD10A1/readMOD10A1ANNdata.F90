!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
! 
! !ROUTINE: readMOD10A1ANNdata
! \label{readMOD10A1ANNdata}
! 
! !REVISION HISTORY: 
!  21 July 2010: Sujay Kumar, Initial Specification
! 
! !INTERFACE: 
subroutine readMOD10A1ANNdata(n,iomode,sindex,eindex)
! !USES:   
  use ESMF
  use netcdf
  use LDT_coreMod
  use LDT_timeMgrMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_ANNMod
  use MOD10A1_ANNdataMod
  use map_utils

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: iomode
  integer, intent(in) :: sindex
  integer, intent(in) :: eindex
! 
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for the synthetic
! soil moisture retrieval product. 
!
!EOP

  integer, parameter           :: nc = 36000, nr=15000
  logical                      :: alarmCheck
  logical                      :: file_exists
  integer                      :: c,r,i,j,c1,r1
  character(len=LDT_CONST_PATH_LEN) :: fname
  logical*1, allocatable       :: lb(:)
  character*1                  :: mod10a1(nc,nr)
  real, allocatable            :: snfrac1(:,:)
  real                         :: snfrac(LDT_rc%lnc(n), LDT_rc%lnr(n))
  real                         :: lat, lon
  integer                      :: ftn,ier,ivar1,scfid
  character*3                  :: fnest

  write(fnest,'(i3.3)') n
  alarmCheck = LDT_isAlarmRinging(LDT_rc,"MOD10A1 data alarm "//trim(fnest))

  if(alarmCheck) then

     call create_ANNMOD10A1_filename(MOD10A1obs(n)%odir, &
          LDT_rc%yr, LDT_rc%mo, LDT_rc%da, fname)
     
     inquire(file=trim(fname),exist=file_exists)
     if(file_exists) then
        
        allocate(lb(MOD10A1obs(n)%nc*MOD10A1obs(n)%nr))
        allocate(snfrac1(MOD10A1obs(n)%nc,MOD10A1obs(n)%nr))


        write(LDT_logunit,*) '[INFO] Reading ..',trim(fname)
        
        ier = nf90_open(path=trim(fname), mode=NF90_WRITE,&
             ncid = ftn)
        if(ier.eq.0) then 
           ivar1 = nf90_inq_varid(ftn,'scf',scfid)
           if(ivar1.eq.0) then 
              lat = MOD10A1obs(n)%gridDesc(4)
              lon = MOD10A1obs(n)%gridDesc(5)
              r1 = nint((lat+59.995)/0.01) + 1
              c1 = nint((lon+179.995)/0.01)+1
              
              call LDT_verify(nf90_get_var(ftn,scfid,snfrac1,&
                   start=(/c1,r1/), &
                   count=(/MOD10A1obs(n)%nc, MOD10A1obs(n)%nr/)),&
                   'Error in nf90_get_var: scf')
           else
              snfrac1 = LDT_rc%udef
           endif
        else
           snfrac1 = LDT_rc%udef
        endif
        call LDT_verify(nf90_close(ftn), 'Error in nf90_close')

        lb = .true. 

        do r=1,MOD10A1obs(n)%nr
           do c=1,MOD10A1obs(n)%nc
              if(snfrac1(c,r).lt.0.or.snfrac1(c,r).gt.100)  then 
                 lb(c+(r-1)*MOD10A1obs(n)%nc) = .false. 
              endif
           enddo
        enddo
       
        call interp_mod10a1(n,MOD10A1obs(n)%nc, &
             MOD10A1obs(n)%nr, &
             snfrac1, lb, snfrac)

        deallocate(lb)
        deallocate(snfrac1)
        
     else
        snfrac = -9999.0
     endif
  else
     snfrac = -9999.0
  endif
     
  call LDT_logSingleANNdata(n,&
       snfrac,  &
       pindex=sindex, &
       iomode = iomode, &
       name = "Snowcover",&
       units="-")
  
end subroutine readMOD10A1ANNdata


!BOP
! 
! !ROUTINE: interp_mod10a1
!  \label{interp_mod10a1}
!
! !INTERFACE: 
subroutine interp_mod10a1(n, nc, nr, var_input, lb, var_output)
! 
! !USES:   
  use LDT_coreMod
  use MOD10A1_ANNdataMod, only : mod10a1obs

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
  integer            :: n
  integer            :: nc
  integer            :: nr
  real               :: var_input(nc*nr)
  logical*1          :: lb(nc*nr)
  real               :: var_output(LDT_rc%lnc(n), LDT_rc%lnr(n))
!EOP
  integer            :: iret
  integer            :: c,r
  logical*1          :: lo(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real               :: go(LDT_rc%lnc(n)*LDT_rc%lnr(n))

  var_output = LDT_rc%udef
  call upscaleByAveraging(&
       mod10a1obs(n)%nc*mod10a1obs(n)%nr, &
       LDT_rc%lnc(n)*LDT_rc%lnr(n), LDT_rc%udef, &
       mod10a1obs(n)%n11, lb, &
       var_input, lo, go)
  
  do r=1, LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
        if(go(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0) then 
           var_output(c,r) = go(c+(r-1)*LDT_rc%lnc(n))/100.0
        endif
     enddo
  enddo
  
end subroutine interp_mod10a1

!BOP
! !ROUTINE: create_ANNMOD10A1_filename
! \label{create_ANNMOD10A1_filename}
! 
! !INTERFACE: 
subroutine create_ANNMOD10A1_filename(ndir, yr, mo,da, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the synthetic filename based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the synthetic soil moisture directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated synthetic filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
 
  filename = trim(ndir)//'/'//&
       trim(fyr)//'/MOD10A1_'//trim(fyr)//trim(fmo)//trim(fda)//'_c5_1km.nc4'
  
end subroutine create_ANNMOD10A1_filename
