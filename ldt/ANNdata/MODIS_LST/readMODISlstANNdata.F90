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
! !ROUTINE: readMODISlstANNdata
! \label{readMODISlstANNdata}
! 
! !REVISION HISTORY: 
!  21 July 2010: Sujay Kumar, Initial Specification
! 
! !INTERFACE: 
subroutine readMODISlstANNdata(n,iomode,sindex,eindex)
! !USES:   
  use ESMF
  use netcdf
  use LDT_coreMod
  use LDT_timeMgrMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_ANNMod
  use MODISlst_ANNdataMod
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

  real              :: timenow
  logical           :: alarmCheck
  logical           :: file_exists
  integer           :: c,r,i,j,c1,r1
  integer           :: ftn
  character(len=LDT_CONST_PATH_LEN) :: fname
  character*3       :: fnest
  integer           :: ier,ivar1, ivar2, ivar3
  integer           :: lstid, mainqcid, lsterrid
  real              :: lat1, lon1
  real              :: lst(MODISlstobs(n)%nc,MODISlstobs(n)%nr)
  real              :: mainqc(MODISlstobs(n)%nc,MODISlstobs(n)%nr)
  real              :: lsterr(MODISlstobs(n)%nc,MODISlstobs(n)%nr)
  real              :: lst1d(MODISlstobs(n)%nc*MODISlstobs(n)%nr)
  logical*1         :: li(MODISlstobs(n)%nc*MODISlstobs(n)%nr)
  real              :: lstobs(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  logical*1         :: lo(LDT_rc%lnc(n)*LDT_rc%lnr(n))

  write(fnest,'(i3.3)') n
  alarmCheck = LDT_isAlarmRinging(LDT_rc,"MODIS LST data alarm "//trim(fnest))

  lstobs= LDT_rc%udef        

  if(alarmCheck) then

     call create_ANNMODISlst_filename(MODISlstobs(n)%odir, &
          LDT_rc%yr, LDT_rc%mo, LDT_rc%da, fname)
     
     inquire(file=trim(fname),exist=file_exists)
     if(file_exists) then
        
        write(LDT_logunit,*) '[INFO] Reading ..',trim(fname)
        
        ier = nf90_open(path=trim(fname), mode=NF90_NOWRITE, &
             ncid= ftn)
        if(ier.eq.0) then 
           ivar1 = nf90_inq_varid(ftn,'lst_daytime',lstid)
           ivar2 = nf90_inq_varid(ftn,'mainqc_day',mainqcid)
           ivar3 = nf90_inq_varid(ftn,'lsterr_day',lsterrid)
           
           if(ivar1.eq.0.and.ivar2.eq.0.and.ivar3.eq.0) then 

              lat1 = MODISlstobs(n)%gridDesc(4)
              lon1 = MODISlstobs(n)%gridDesc(5)
              r1 = nint((lat1+59.995)/0.01) + 1
              c1 = nint((lon1+179.995)/0.01) + 1

              call LDT_verify(nf90_get_Var(ftn,lstid, lst,&
                   start=(/c1,r1/), &
                   count=(/MODISlstobs(n)%nc, MODISlstobs(n)%nr/)),&
                   'Error in nf90_get_var: lst_daytime')
              call LDT_verify(nf90_get_Var(ftn,mainqcid, mainqc,&
                   start=(/c1,r1/), &
                   count=(/MODISlstobs(n)%nc, MODISlstobs(n)%nr/)),&
                   'Error in nf90_get_var: mainqc_day')
              call LDT_verify(nf90_get_Var(ftn,lsterrid, lsterr,&
                   start=(/c1,r1/), &
                   count=(/MODISlstobs(n)%nc, MODISlstobs(n)%nr/)),&
                   'Error in nf90_get_var: lsterr_day')
           else
              lst    = LDT_rc%udef
              mainqc = 3
              lsterr = 3
           endif           
        else
           lst    = LDT_rc%udef
           mainqc = 3
           lsterr = 3
        endif
        
        call LDT_verify(nf90_close(ftn), 'Error in nf90_close')

        li = .false. 
        
!values
        do r=1,MODISlstobs(n)%nr
           do c=1,MODISlstobs(n)%nc
              if(lst(c,r).ne.-9999.0.and.&
                   (mainqc(c,r).eq.0.or.mainqc(c,r).eq.1).and.&
                   (lsterr(c,r).eq.0)) then 
                 lst1d(c+(r-1)*MODISlstobs(n)%nc) = lst(c,r) 
                 li(c+(r-1)*MODISlstobs(n)%nc) = .true. 
              else
                 lst1d(c+(r-1)*MODISlstobs(n)%nc) = -9999.0
                 li(c+(r-1)*MODISlstobs(n)%nc) = .false. 
              endif
           enddo
        enddo
        
        call upscaleByAveraging(MODISlstobs(n)%nc*MODISlstobs(n)%nr,&
             LDT_rc%lnc(n)*LDT_rc%lnr(n),LDT_rc%udef,&
             MODISlstobs(n)%n11,li,lst1d,lo,lstobs)

!        do r=1,LDT_rc%lnr(n)
!           do c=1,LDT_rc%lnc(n)
!              if(lstobs(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0) then 
!                 MODISlstobs(n)%lstobs(c,r) = lstobs(c+(r-1)*LDT_rc%lnc(n))
!              endif
!           enddo
!        enddo
        
     endif
  endif

  call LDT_logSingleANNdata(n,&
       lstobs,  &
       pindex=sindex, &
       iomode = iomode, &
       name = "AvgSurfT",&
       units="K")
  
end subroutine readMODISlstANNdata

!BOP
! !ROUTINE: create_ANNMODISlst_filename
! \label{create_ANNMODISlst_filename}
! 
! !INTERFACE: 
subroutine create_ANNMODISlst_filename(ndir, yr, mo,da, filename)
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
       trim(fyr)//'.'//trim(fmo)//'.'//trim(fda)//'/LST_Day_1km_all.nc4'
  
end subroutine create_ANNMODISlst_filename
