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
! !ROUTINE: readMODIS_LSTObs
! \label{readMODIS_LSTObs}
!
! !INTERFACE: 
subroutine readMODIS_LSTObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,     only : LVT_rc
  use LVT_logMod,      only : LVT_logunit, LVT_verify, &
       LVT_getNextUnitNumber, LVT_releaseUnitNumber
  use LVT_histDataMod
  use MODIS_LSTobsMod, only : MODISLSTObs
          
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in)    :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!  NOTES: 
!   The MODIS LST data is available at daily intervals. So 
!   the comparisons against model data should use at least a 
!   24 hour (1day) averaging interval. 
!
!  Main qc flag 
!   0 - LST produced, good quality, not necessary to examine more detailed QA
!   1 - LST produced, other quality, recommend examination of more detailed QA
!   2 - LST not produced due to cloud effects
!   3 - LST not produced primarily due to reasons other than cloud
!   
!  LST error flag key: 
!   0 = average LST error <=1k
!   1 = average LST error <=2k
!   2 = average LST error <=3k
!   3 = average LST error > 3k
!
!  Emissivity error flag key: 
!   0 = average LST error <=0.01
!   1 = average LST error <=0.02
!   2 = average LST error <=0.04
!   3 = average LST error > 0.04
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  7 Mar 2015: Sujay Kumar, Initial Specification
! 
!EOP

  character*100          :: filename
  integer                :: ftn 
  logical                :: file_exists
  integer                :: nid, ios, ier, ivar1, ivar2, ivar3
  real                   :: time
  integer                :: c1,r1,line
  real                   :: lst(modislstobs(source)%nc,modislstobs(source)%nr)
  real                   :: mainqc(modislstobs(source)%nc,modislstobs(source)%nr)
  real                   :: lsterr(modislstobs(source)%nc,modislstobs(source)%nr)
  real                   :: lst1d(modislstobs(source)%nc*modislstobs(source)%nr)
  logical*1              :: li(modislstobs(source)%nc*modislstobs(source)%nr)
  real                   :: lat1,lon1
  integer                :: c,r,t,kk
  integer                :: lstid, mainqcid, lsterrid
  logical*1              :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real                   :: varfield(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: gridDesc(6)

  
  Modislstobs(Source)%lst = LVT_rc%udef
  varfield = LVT_rc%udef

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  
  call create_modislst_filename(Modislstobs(Source)%odir, &
       LVT_rc%dyr(source),LVT_rc%dmo(source),LVT_rc%dda(source), filename)
  
  time = LVT_rc%dhr(source)*3600+LVT_rc%dmn(source)*60+LVT_rc%dss(source)
  if(mod(time,86400.0).eq.0.0.or.LVT_rc%resetFlag(source)) then 

     LVT_rc%resetFlag(source) = .false. 
     inquire(file=trim(filename),exist=file_exists) 
     if(file_exists) then 
        write(LVT_logunit,*) '[INFO] Reading MODIS LST file ',trim(filename)
        
        gridDesc = 0 
        gridDesc(1) = modislstobs(source)%gridDesc(4)
        gridDesc(2) = modislstobs(source)%gridDesc(5)
        gridDesc(3) = modislstobs(source)%gridDesc(7)
        gridDesc(4) = modislstobs(source)%gridDesc(8)
        gridDesc(5) = 0.01
        gridDesc(6) = 0.01
        
        ier = nf90_open(path=trim(filename),mode=NF90_NOWRITE, &
             ncid=ftn)
        
        if(ier.eq.0) then 
           ivar1 = nf90_inq_varid(ftn,'lst_daytime',lstid)
           ivar2 = nf90_inq_varid(ftn,'mainqc_day',mainqcid)
           ivar3 = nf90_inq_varid(ftn,'lsterr_day',lsterrid)
           
           if(ivar1.eq.0.and.ivar2.eq.0.and.ivar3.eq.0) then 
              lat1 = modislstobs(source)%gridDesc(4)
              lon1 = modislstobs(source)%gridDesc(5)
              r1 = nint((lat1+59.995)/0.01)+1
              c1 = nint((lon1+179.995)/0.01)+1
              
              call LVT_verify(nf90_get_var(ftn,lstid,lst, &
                   start=(/c1,r1/),&
                   count=(/modislstobs(source)%nc, modislstobs(source)%nr/)),&
                   'Error in nf90_get_var: lst_daytime')
              
              call LVT_verify(nf90_get_var(ftn,mainqcid,mainqc, &
                   start=(/c1,r1/),&
                   count=(/modislstobs(source)%nc, modislstobs(source)%nr/)),&
                   'Error in nf90_get_var: mainqc_day')
              
              call LVT_verify(nf90_get_var(ftn,lsterrid,lsterr, &
                   start=(/c1,r1/),&
                   count=(/modislstobs(source)%nc, modislstobs(source)%nr/)),&
                   'Error in nf90_get_var: lsterr_day')
           else
              lst    = LVT_rc%udef
              mainqc = 3
              lsterr = 3
           endif
        else
           lst    = LVT_rc%udef
           mainqc = 3
           lsterr = 3
        endif
        call LVT_verify(nf90_close(ftn),&
             'Error in nf90_close')
        
        li = .false. 
        !values
        do r=1,modislstobs(source)%nr
           do c=1,modislstobs(source)%nc
              if(lst(c,r).ne.-9999.0.and.&
                   (mainqc(c,r).eq.0.or.mainqc(c,r).eq.1).and.&
                   (lsterr(c,r).eq.0)) then 
                 lst1d(c+(r-1)*modislstobs(source)%nc) = lst(c,r) 
                 li(c+(r-1)*modislstobs(source)%nc) = .true. 
              else
                 lst1d(c+(r-1)*modislstobs(source)%nc) = -9999.0
                 li(c+(r-1)*modislstobs(source)%nc) = .false. 
              endif
           enddo
        enddo
        
        call upscaleByAveraging(modislstobs(source)%nc*modislstobs(source)%nr,&
             LVT_rc%lnc*LVT_rc%lnr,LVT_rc%udef,&
             Modislstobs(Source)%n11,li,lst1d,lo,Modislstobs(Source)%lst)
        
        do r=1, LVT_rc%lnr
           do c=1, LVT_rc%lnc
              if(lo(c+(r-1)*LVT_rc%lnc)) then 
                 varfield(c,r) = Modislstobs(source)%lst(c+(r-1)*LVT_rc%lnc) 
              else
                 varfield(c,r) = LVT_rc%udef
              endif
           enddo
        enddo
        
     else
        varfield  = LVT_rc%udef
     endif
  endif
#endif
  
  call LVT_logSingleDataStreamVar(LVT_MOC_AVGSURFT,source,varfield,vlevel=1,units="K")
  call LVT_logSingleDataStreamVar(LVT_MOC_RADT,source,varfield,vlevel=1,units="K")
  
end subroutine readMODIS_LSTObs

!BOP
! 
! !ROUTINE: create_modislst_filename
! \label{create_modislst_filename}
!
! !INTERFACE: 
subroutine create_modislst_filename(odir,yr,mo,da,filename)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for MODIS LST data files 
! based on the given date (year, month, day)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      GlboSnow base directory
!   \item[yr]        year of data
!   \item[filename]  Name of the MODISLST_LH file
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  character(len=*)             :: odir
  integer                      :: yr
  integer                      :: mo
  integer                      :: da
  character(len=*)             :: filename
!
!EOP

  character*4             :: fyr
  character*2             :: fmo
  character*2             :: fda

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  
  filename = trim(odir)//'/'//trim(fyr)//'.'//trim(fmo)//'.'//&
       trim(fda)//'/LST_Day_1km_all.nc4'
  
end subroutine create_modislst_filename


