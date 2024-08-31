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
! !ROUTINE: readLISDAdiagOutput
! \label(readLISDAdiagOutput)
!
! !INTERFACE:
subroutine readLISDAdiagOutput(source)
! !USES:   
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif
  use LVT_coreMod,      only : LVT_rc
  use LVT_histDataMod
  use LVT_logMod
  use LISDAdiagOutputMod,    only : lisdadiagoutput

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in) :: source 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!
! !NOTES: Currently limited to the LIS6 outputs produced for the 
! NOHRSC with Noah 3.2
! 
!EOP

!  character*100     :: fname 
  character*140     :: fname 
  logical           :: file_exists
  real              :: spread_var(lisdadiagoutput(source)%nstvars,&
       LVT_rc%lnc,LVT_rc%lnr)
  real              :: incr_var(lisdadiagoutput(source)%nstvars,&
       LVT_rc%lnc,LVT_rc%lnr)
  real              :: innov_var(LVT_rc%lnc,LVT_rc%lnr)
  real              :: obscount_var(LVT_rc%lnc,LVT_rc%lnr)
  integer           :: varid
  integer           :: t,iret
  integer           :: ftn
  integer           :: c,r

  character*2        :: cda
  character*100      :: cdate, cdate1
  character(len=100) :: var_name
  integer            :: k
   
  spread_var = LVT_rc%udef
  incr_var  = LVT_rc%udef

  write(unit=cda, fmt='(i2.2)') &
       lisdadiagoutput(source)%instance
  
  if(lisdadiagoutput(source)%computeSpread.eq.1) then 
     write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
          LVT_rc%dyr(source), LVT_rc%dmo(source), LVT_rc%dda(source), &
          LVT_rc%dhr(source), LVT_rc%dmn(source)
     
     write(unit=cdate, fmt='(i4.4, i2.2)') &
          LVT_rc%dyr(source), LVT_rc%dmo(source)
     
     fname = trim(lisdadiagoutput(source)%odir)//'/EnKF/'//trim(cdate)//'/'&
          //'LIS_DA_EnKF_'//trim(cdate1)//'_spread.a'//trim(cda)//'.d01.nc'
     
     inquire(file=trim(fname),exist=file_exists)
     
     if(file_exists) then
        write(LVT_logunit,*) '[INFO] reading LIS DA ensemble spread output ',&
             trim(fname)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        iret = nf90_open(path=trim(fname),mode=NF90_NOWRITE, &
             ncid = ftn)
        if(iret.eq.0) then 
           
           do k=1,lisdadiagoutput(source)%nstvars
              iret = nf90_inquire_variable(ftn, k, var_name)
              call LVT_verify(nf90_inq_varid(ftn,trim(var_name),varid),&
                   'nf90_inq_varid failed for '//trim(var_name))
              call LVT_verify(nf90_get_var(ftn,varid,spread_var(k,:,:)),&
                   'Error in nf90_get_var for '//trim(var_name))
           enddo
        endif
        iret = nf90_close(ftn)
     else
        write(LVT_logunit,*) '[WARN] Warning: LIS DA ensemble spread file ',&
             trim(fname),' does not exist'
        spread_var = -9999.0
     endif
#endif
     
     do k=1,lisdadiagoutput(source)%nstvars
        call LVT_logSingleDataStreamVar(LVT_MOC_DA_ENSSPREAD, source,&
             spread_var(k,:,:),vlevel=k,units="-")
     enddo
  endif

  if(lisdadiagoutput(source)%computeAnlIncr.eq.1) then 
     write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
          LVT_rc%dyr(source), LVT_rc%dmo(source), &
          LVT_rc%dda(source), LVT_rc%dhr(source), LVT_rc%dmn(source)
     
     write(unit=cdate, fmt='(i4.4, i2.2)') &
          LVT_rc%dyr(source), LVT_rc%dmo(source)
     
     fname = trim(lisdadiagoutput(source)%odir)//'/EnKF/'//trim(cdate)//'/'&
          //'LIS_DA_EnKF_'//trim(cdate1)//'_incr.a'//trim(cda)//'.d01.nc'
     
     inquire(file=trim(fname),exist=file_exists)
     
     if(file_exists) then
        write(LVT_logunit,*) '[INFO] reading LIS DA ensemble increments output ',&
             trim(fname)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        iret = nf90_open(path=trim(fname),mode=NF90_NOWRITE, &
             ncid = ftn)
        if(iret.eq.0) then 
           
           do k=1,lisdadiagoutput(source)%nstvars
              iret = nf90_inquire_variable(ftn, k, var_name)
              call LVT_verify(nf90_inq_varid(ftn,trim(var_name),varid),&
                   'nf90_inq_varid failed for '//trim(var_name))
              call LVT_verify(nf90_get_var(ftn,varid, incr_var(k,:,:)),&
                   'Error in nf90_get_var for '//trim(var_name))
           enddo
        endif
        iret = nf90_close(ftn)

     else
        write(LVT_logunit,*) '[WARN] Warning: LIS DA ensemble increments file ',&
             trim(fname),' does not exist'
        incr_var = -9999.0
     endif
#endif
     
     do k=1,lisdadiagoutput(source)%nstvars
        call LVT_logSingleDataStreamVar(LVT_MOC_DA_INCR, source,&
             incr_var(k,:,:),vlevel=k,units="-")
     enddo
  endif

  innov_var = LVT_rc%udef
  obscount_var = LVT_rc%udef

  if(lisdadiagoutput(source)%computeInnovDist.eq.1) then 
     write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
          LVT_rc%dyr(source), LVT_rc%dmo(source), LVT_rc%dda(source), &
          LVT_rc%dhr(source), LVT_rc%dmn(source)
     
     write(unit=cdate, fmt='(i4.4, i2.2)') &
          LVT_rc%dyr(source), LVT_rc%dmo(source)
     
     fname = trim(lisdadiagoutput(source)%odir)//'/EnKF/'//trim(cdate)//'/'&
          //'LIS_DA_EnKF_'//trim(cdate1)//'_innov.a'//trim(cda)//'.d01.nc'
     
     inquire(file=trim(fname),exist=file_exists)
     
     if(file_exists) then
        write(LVT_logunit,*) '[INFO] reading LIS DA ensemble innov output ',&
             trim(fname)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        iret = nf90_open(path=trim(fname),mode=NF90_NOWRITE, &
             ncid = ftn)
        if(iret.eq.0) then 
           
           var_name = 'ninnov_'//trim(cda)
           call LVT_verify(nf90_inq_varid(ftn,trim(var_name),varid),&
                'nf90_inq_varid failed for '//trim(var_name))
           call LVT_verify(nf90_get_var(ftn,varid,innov_var(:,:)),&
                'Error in nf90_get_var for '//trim(var_name))

           do r=1,LVT_rc%lnr
              do c=1,LVT_rc%lnc
                 if(innov_var(c,r).ne.LVT_rc%udef) then 
                    obscount_var(c,r) = 1
                 endif
              enddo
           enddo

        endif
        iret = nf90_close(ftn)
     else
        write(LVT_logunit,*) '[WARN] Warning: LIS DA ensemble innov file ',&
             trim(fname),' does not exist'
        innov_var = -9999.0
     endif
#endif

     call LVT_logSingleDataStreamVar( LVT_MOC_DA_NINNOV, source,&
          innov_var(:,:),vlevel=1,units="-")
     call LVT_logSingleDataStreamVar( LVT_MOC_DA_OBSCOUNT, source,&
          obscount_var(:,:),vlevel=1,units="-")
  endif

end subroutine readLISDAdiagOutput

